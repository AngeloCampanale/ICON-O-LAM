! Contains the implementation of the boundary conditions for the limited area mode
! for the ocean (ICON-O-LAM).
!
! This module controls the boundary conditions as well as ...?
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_ocean_limarea
  !-------------------------------------------------------------------------
  ! USE section
  !-------------------------------------------------------------------------
  USE mo_kind,              ONLY: wp, i8
  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_model_domain,      ONLY: t_patch, t_patch_3d, p_patch
  USE mo_sync,              ONLY: sync_c, sync_e, sync_patch_array
  USE mo_read_interface,    ONLY: read_2D_1time, read_3D_1time, on_cells, on_edges, t_stream_id, &
    & openInputFile, closeFile
  USE mo_grid_subset,       ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,   ONLY: nproma
  USE mo_ocean_nml,         ONLY: n_zlev, ocean_latbc_path, ocean_latbc_dtime
  USE mo_ocean_nml,         ONLY: ocean_latbc_filepattern, ocean_latbc_from_vn
  USE mo_io_units,          ONLY: filename_max
  USE mo_impl_constants,    ONLY: SUCCESS
  USE mtime,                ONLY: timedelta, datetime, getPTStringFromMS, &
    &                             newTimedelta, MAX_TIMEDELTA_STR_LEN, &
    &                             MAX_DATETIME_STR_LEN, getTotalSecondsTimeDelta, &
    &                             OPERATOR(+), OPERATOR(-), OPERATOR(>)
  USE mo_util_mtime,        ONLY: mtime_utils, FMT_DDDHH, FMT_DDHHMMSS, FMT_HHH
  USE mo_impl_constants,    ONLY: MAX_CHAR_LENGTH
  USE mo_util_string,       ONLY: t_keyword_list, int2string, &
    &                             associate_keyword, with_keywords
  USE mo_grid_config,       ONLY: nroot, n_dom
  USE mo_time_config,       ONLY: time_config
  USE mo_mpi,               ONLY: get_my_global_mpi_id
  USE mo_ocean_types,       ONLY: t_hydro_ocean_state
  USE mo_dynamics_config,   ONLY: nold, nnew
  USE mo_math_types,        ONLY: t_cartesian_coordinates
  USE mo_scalar_product,    ONLY: map_cell2edges_3D
  USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
  USE mo_math_utilities,    ONLY: gvec2cvec
  USE mo_math_constants,    ONLY: rad2deg! , pi, pi_2, deg2rad
!  USE mo_fortran_tools,     ONLY: assign_if_present

  IMPLICIT NONE
  PRIVATE

  !-------------------------------------------------------------------------
  ! PUBLIC declarations
  !-------------------------------------------------------------------------
  ! Attributes
  PUBLIC :: ocean_latbc_data

  ! Methods
  PUBLIC :: init_ocean_latbc
  PUBLIC :: destruct_ocean_latbc
  PUBLIC :: update_ocean_latbc
  PUBLIC :: apply_ocean_ssh_latbc
  PUBLIC :: apply_ocean_velocity_latbc
  PUBLIC :: apply_ocean_tracer_latbc

  !-------------------------------------------------------------------------
  ! TYPE definitions
  !-------------------------------------------------------------------------
  TYPE t_ocean_latbc_data ! make sure to test if stretch_c is needed
    REAL(wp), ALLOCATABLE :: to(:,:,:,:)
    REAL(wp), ALLOCATABLE :: so(:,:,:,:)
    REAL(wp), ALLOCATABLE :: eta_c(:,:,:)
    REAL(wp), ALLOCATABLE :: stretch_c(:,:,:)
    REAL(wp), ALLOCATABLE :: vn(:,:,:,:)
    REAL(wp), ALLOCATABLE :: u(:,:,:,:)
    REAL(wp), ALLOCATABLE :: v(:,:,:,:)

    REAL(wp), ALLOCATABLE :: nudge_cell(:,:)
    REAL(wp), ALLOCATABLE :: nudge_edge(:,:)
  CONTAINS
    PROCEDURE :: alloc => t_ocean_latbc_data_allocate
    PROCEDURE :: finalize => t_ocean_latbc_data_finalize
  END TYPE t_ocean_latbc_data

  !-------------------------------------------------------------------------
  ! MODULE PARAMETERS
  !-------------------------------------------------------------------------
  CHARACTER(LEN=16), PARAMETER :: module_name = 'mo_ocean_limarea'

  INTEGER, PARAMETER :: OCEAN_LATBC_TYPE_CONST       = 0
  INTEGER, PARAMETER :: OCEAN_LATBC_TYPE_VAR         = 1

  !-------------------------------------------------------------------------
  ! MODULE VARIABLES
  !-------------------------------------------------------------------------
  TYPE(t_ocean_latbc_data), TARGET :: ocean_latbc_data
  INTEGER :: ocean_latbc_timelevel ! This is the timelevel of younger latbc. We save only this index, the other is calculated on the fly

  TYPE(timedelta), POINTER            :: ocean_latbc_dtime_mtime ! dt between latbc levels
  CHARACTER(LEN=filename_max)         :: ocean_latbc_filename ! filename read by the module methods
  TYPE(datetime)                      :: ocean_latbc_prev_datetime ! mtime date of older latbc data
  TYPE(datetime)                      :: ocean_latbc_next_datetime ! mtime date of newer latbc data
  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: ocean_latbc_next_datetime_str ! string version

  ! TYPE(t_operator_coeff), POINTER :: this_operators_coeff
  ! TYPE(t_hydro_ocean_state), POINTER :: this_ocean_state

  ! CONSTANTS (for better readability):
  ! OCEAN_LATBC_TYPE_CONST: constant lateral boundary conditions
  ! OCEAN_LATBC_TYPE_VAR: time-dependent lateral boundary conditions

CONTAINS

  SUBROUTINE init_ocean_latbc(patch_3d, operators_coeff)
    TYPE(t_patch_3d), POINTER, INTENT(in) :: patch_3d
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff

    TYPE(t_subset_range), POINTER :: all_cells, all_edges
    TYPE(t_patch), POINTER :: patch_2d
    INTEGER :: numCellBlocks, numEdgeBlocks
    REAL(wp) :: ocean_latbc_dtime_in_ms
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: ocean_latbc_dtime_str
    INTEGER :: errno
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_index_e, end_index_e
    INTEGER :: jc, je, blockNo
    CHARACTER(*), PARAMETER :: method_name = module_name//"::init_ocean_latbc"

    patch_2d => patch_3d%p_patch_2d(1)

    ! allocate the arrays for latbc data
    numCellBlocks = patch_2d%alloc_cell_blocks
    numEdgeBlocks = patch_2d%nblks_e
    CALL ocean_latbc_data%alloc(numCellBlocks, numEdgeBlocks)

    ! convert ocean_latbc_dtime (integer) into mtime object (datetime)
    ocean_latbc_dtime_in_ms = 1000._wp * ocean_latbc_dtime
    CALL getPTStringFromMS(NINT(ocean_latbc_dtime_in_ms,i8), ocean_latbc_dtime_str)
    ocean_latbc_dtime_mtime => newTimedelta(ocean_latbc_dtime_str, errno)
    IF (errno /= 0)  CALL finish(method_name, "Error in initialization of dtime_latbc time delta.")

    ! prepare some debugging
    !CALL save_const_zone_data(patch_2d)

    ! now preload the first two levels
    ! set the datetime and filename
    ocean_latbc_next_datetime = time_config%tc_exp_startdate ! does this always work?
    ocean_latbc_filename = TRIM(ocean_latbc_path) &
      & // TRIM(generate_filename_ocean(ocean_latbc_next_datetime))
    IF (get_my_global_mpi_id() == 0) write(*,*) "reading 1st latbc file: ", TRIM(ocean_latbc_filename)
    CALL read_3D_latbc_from_file(patch_3d, 1, operators_coeff)
    !CALL save_global(0, 1, patch_2d)

    ! set the datetimes, filename and also the timelevel
    ocean_latbc_prev_datetime = ocean_latbc_next_datetime
    ocean_latbc_next_datetime = ocean_latbc_next_datetime + ocean_latbc_dtime_mtime
    ocean_latbc_filename = TRIM(ocean_latbc_path) &
      & // TRIM(generate_filename_ocean(ocean_latbc_next_datetime))
    IF (get_my_global_mpi_id() == 0) write(*,*) "reading 2nd latbc file: ", TRIM(ocean_latbc_filename)
    CALL read_3D_latbc_from_file(patch_3d, 2, operators_coeff)
    !CALL save_global(1, 2, patch_2d)
    ocean_latbc_timelevel = 2

    ! prepare nudging coefficients
    ! nudge cells
    all_cells => patch_2d%cells%all
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c ! do it later using ordering
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 4 .and. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= 8) THEN ! later do it properly, using namelist parameters (also crosscheck the used doesnt require too much)
          ocean_latbc_data%nudge_cell(jc,blockNo) = &
            & 0.5*exp(-(patch_2d%cells%refin_ctrl(jc,blockNo)-4)/2._wp)
        END IF
      END DO ! jc
    END DO ! blockNo

    ! nudge edges
    all_edges => patch_2d%edges%all
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
      DO je = start_index_e, end_index_e ! do it later using oredering
        IF (patch_2d%edges%refin_ctrl(je,blockNo) > 8 .and. &
          & patch_2d%edges%refin_ctrl(je,blockNo) <= 16) THEN ! later do it properly, using namelist parameters
          ocean_latbc_data%nudge_edge(je,blockNo) = &
            & 0.5*exp(-(patch_2d%edges%refin_ctrl(je,blockNo)-8)/4._wp)
        END IF
      END DO ! je
    END DO ! blockNo

    ! mapping for sparse grid would also be done here if it was there

  END SUBROUTINE init_ocean_latbc

  SUBROUTINE destruct_ocean_latbc()
    CALL ocean_latbc_data%finalize
  END SUBROUTINE destruct_ocean_latbc

  FUNCTION generate_filename_ocean(latbc_mtime, opt_mtime_begin) RESULT(result_str)
    CHARACTER(MAX_CHAR_LENGTH)                       :: result_str
    TYPE(datetime),   INTENT(IN)                     :: latbc_mtime
    ! Optional: Start date, which a time span in the filename is related to.
    TYPE(datetime),   INTENT(IN),  POINTER, OPTIONAL :: opt_mtime_begin
    ! Local variables
    CHARACTER(MAX_CHAR_LENGTH), PARAMETER :: method_name = module_name//'::generate_filename_ocean'
    TYPE (t_keyword_list), POINTER        :: keywords => NULL()
    CHARACTER(MAX_CHAR_LENGTH)            :: str
    INTEGER :: jlev

    jlev = p_patch(1)%level ! not tested yet - does it work in the ocean the same as in the atmo

    WRITE(str,'(i4)')   latbc_mtime%date%year
    CALL associate_keyword("<y>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%month
    CALL associate_keyword("<m>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%date%day
    CALL associate_keyword("<d>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%hour
    CALL associate_keyword("<h>",         TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%minute
    CALL associate_keyword("<min>",       TRIM(str),                        keywords)
    WRITE(str,'(i2.2)') latbc_mtime%time%second !FLOOR(latbc_mtime%time%second)
    CALL associate_keyword("<sec>",       TRIM(str),                        keywords)

    CALL associate_keyword("<nroot>",     TRIM(int2string(nroot,'(i1)')),   keywords)
    CALL associate_keyword("<nroot0>",    TRIM(int2string(nroot,'(i2.2)')), keywords)
    CALL associate_keyword("<jlev>",      TRIM(int2string(jlev, '(i2.2)')), keywords)
    CALL associate_keyword("<dom>",       TRIM(int2string(1,'(i2.2)')),     keywords)

    IF (PRESENT(opt_mtime_begin)) THEN
      CALL associate_keyword("<ddhhmmss>", &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDHHMMSS)), &
        &                    keywords)
      CALL associate_keyword("<dddhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_DDDHH)),    &
        &                    keywords)
      CALL associate_keyword("<hhh>",    &
        &                    TRIM(mtime_utils%ddhhmmss(opt_mtime_begin, latbc_mtime, FMT_HHH)),    &
        &                    keywords)
    END IF

    ! replace keywords in latbc_filename
    result_str = TRIM(with_keywords(keywords, TRIM(ocean_latbc_filepattern)))
  END FUNCTION generate_filename_ocean

  SUBROUTINE t_ocean_latbc_data_allocate(this, numCellBlocks, numEdgeBlocks)
    CLASS(t_ocean_latbc_data) :: this
    INTEGER, INTENT(in) :: numCellBlocks
    INTEGER, INTENT(in) :: numEdgeBlocks
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//'::allocate_ocean_latbc'
    INTEGER :: ierrstat=0

    ALLOCATE(this%to(3,nproma,n_zlev,numCellBlocks))
    ALLOCATE(this%so(3,nproma,n_zlev,numCellBlocks))
    ALLOCATE(this%eta_c(3,nproma,numCellBlocks))
    ALLOCATE(this%stretch_c(3,nproma,numCellBlocks))
    ALLOCATE(this%u(3,nproma,n_zlev,numCellBlocks))
    ALLOCATE(this%v(3,nproma,n_zlev,numCellBlocks))
    ALLOCATE(this%vn(3,nproma,n_zlev,numEdgeBlocks))

    ALLOCATE(this%nudge_cell(nproma,numCellBlocks))
    ALLOCATE(this%nudge_edge(nproma,numEdgeBlocks), STAT=ierrstat)

    IF (ierrstat /= SUCCESS) CALL finish(method_name, "ALLOCATE failed!")
  END SUBROUTINE t_ocean_latbc_data_allocate

  SUBROUTINE t_ocean_latbc_data_finalize(this)
    CLASS(t_ocean_latbc_data) :: this
    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//'::t_ocean_latbc_data_finalize'
    INTEGER :: ierrstat=0

    !CALL message(method_name, 't_ocean_latbc_data_finalize')
    IF (ALLOCATED(this%to)) DEALLOCATE(this%to, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%so)) DEALLOCATE(this%so, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%eta_c)) DEALLOCATE(this%eta_c, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%stretch_c)) DEALLOCATE(this%stretch_c, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%u)) DEALLOCATE(this%u, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%v)) DEALLOCATE(this%v, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%vn)) DEALLOCATE(this%vn, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")

    IF (ALLOCATED(this%nudge_cell)) DEALLOCATE(this%nudge_cell, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
    IF (ALLOCATED(this%nudge_edge)) DEALLOCATE(this%nudge_edge, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(method_name, "DEALLOCATE failed!")
  END SUBROUTINE t_ocean_latbc_data_finalize

  SUBROUTINE configure_ocean_latbc()
  !--------------------------------------------------------------------------------------
  !  Set up parameters
  !--------------------------------------------------------------------------------------
    CHARACTER(*), PARAMETER :: method_name = module_name//&
      "mo_limarea_config::configure_ocean_latbc"

!!!    !----------------------------------------------------
!!!    ! Sanity check and Prints
!!!    !----------------------------------------------------
!!!
!!!    IF (ocean_latbc_itype == OCEAN_LATBC_TYPE_CONST) THEN
!!!
!!!       WRITE(message_text,'(a)')'Lateral boundary nudging using constant boundary data from the initial conditions.'
!!!       CALL message(TRIM(method_name),message_text)
!!!
!!!    ELSE IF (ocean_latbc_itype == OCEAN_LATBC_TYPE_VAR) THEN
!!!
!!!       WRITE(message_text,'(a)')'Lateral boundary condition using interpolated boundary data.'
!!!       CALL message(TRIM(method_name),message_text)
!!!
!!!!       IF (num_prefetch_proc == 0) THEN
!!!!         WRITE(message_text,'(a)') 'Synchronous latBC input has been disabled'
!!!!         CALL finish(TRIM(method_name),message_text)
!!!!       END IF
!!!
!!!    ELSE
!!!
!!!       WRITE(message_text,'(a,i8)') 'Invalid lateral boundary condition mode:', ocean_latbc_itype
!!!       CALL finish(TRIM(method_name),message_text)
!!!
!!!    END IF
!!!
!!!    ! Check whether an mapping file is provided for prefetching boundary data
!!!    ! calls a finish either when the flag is absent
!!!    !
!!!!    IF (ocean_latbc_config%latbc_varnames_map_file == ' ') THEN
!!!!       WRITE(message_text,'(a)') 'no latbc_varnames_map_file provided.'
!!!!       CALL message(TRIM(method_name),message_text)
!!!!    ENDIF


  END SUBROUTINE configure_ocean_latbc
  !--------------------------------------------------------------------------------------

  SUBROUTINE update_ocean_latbc(patch_3d, ocean_state, current_time, previous_time, operators_coeff, jstep) !next_timestep_to_load)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), INTENT(in) :: current_time ! make it pointer?
    TYPE(datetime), INTENT(in) :: previous_time
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff
    INTEGER, INTENT(in) :: jstep

    TYPE(timedelta), POINTER :: ocean_time_step
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(timedelta) :: current_ocean_latbc_dtime_mtime
    REAL(wp) :: current_ocean_latbc_dtime
    REAL(wp) :: latbc_f ! interpolation factor
    INTEGER :: jg, l1, l2

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)
    ocean_time_step => time_config%tc_dt_model

    ! current_time in code stands for the timepoint towards which we are moving, i.e. n+1

    ! However,  we want the timepoint n, so we use (current time - model timestep).
    ! This logic goes for the update at the beginning of the timestep but not for the refresh
    ! after particular solvers where we do want timepoint n+1.

    DO WHILE (previous_time > ocean_latbc_next_datetime)
      ! set the datetimes, filename and the timelevel
      ocean_latbc_prev_datetime = ocean_latbc_next_datetime
      ocean_latbc_next_datetime = ocean_latbc_next_datetime + ocean_latbc_dtime_mtime
      ocean_latbc_filename = TRIM(ocean_latbc_path) &
        & // TRIM(generate_filename_ocean(ocean_latbc_next_datetime))
      IF (get_my_global_mpi_id() == 0) write(*,*) "reading next latbc file: ", TRIM(ocean_latbc_filename)
      ocean_latbc_timelevel = 3 - ocean_latbc_timelevel ! jump to old level, it's the new one now
      CALL read_3D_latbc_from_file(patch_3d, ocean_latbc_timelevel, operators_coeff) ! overwrite old
      CALL save_global(jstep, ocean_latbc_timelevel, patch_2d)
    END DO

    ! calculate the interpolation coefficient
    current_ocean_latbc_dtime_mtime = previous_time - ocean_latbc_prev_datetime
    current_ocean_latbc_dtime = REAL(getTotalSecondsTimeDelta(current_ocean_latbc_dtime_mtime, &
      &                                                       ocean_latbc_prev_datetime))
    latbc_f = current_ocean_latbc_dtime / ocean_latbc_dtime

    IF (get_my_global_mpi_id() == 0) write(*,*) "latbc_f = ", latbc_f

    ! interpolate
    l1 = 3 - ocean_latbc_timelevel
    l2 = ocean_latbc_timelevel
    ocean_latbc_data%to(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%to(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%to(l2,:,:,:) ! 1 for temperature
    ocean_latbc_data%so(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%so(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%so(l2,:,:,:) ! 2 for salinity
    ocean_latbc_data%eta_c(3,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%eta_c(l1,:,:) + &
    &    latbc_f *ocean_latbc_data%eta_c(l2,:,:)
    ocean_latbc_data%vn(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%vn(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%vn(l2,:,:,:)

    CALL apply_ocean_ssh_latbc(patch_3d, ocean_state, nold(jg))
    CALL apply_ocean_velocity_latbc(patch_3d, ocean_state, nold(jg))
    CALL apply_ocean_tracer_latbc(patch_3d, ocean_state, nold(jg))

    ! advance the time pointer if needed, because we will never need timepoint n
    DO WHILE (current_time > ocean_latbc_next_datetime)
      ! set the datetimes, filename and the timelevel
      ocean_latbc_prev_datetime = ocean_latbc_next_datetime
      ocean_latbc_next_datetime = ocean_latbc_next_datetime + ocean_latbc_dtime_mtime
      ocean_latbc_filename = TRIM(ocean_latbc_path) &
        & // TRIM(generate_filename_ocean(ocean_latbc_next_datetime))
      IF (get_my_global_mpi_id() == 0) write(*,*) "reading next latbc file: ", TRIM(ocean_latbc_filename)
      ocean_latbc_timelevel = 3 - ocean_latbc_timelevel ! jump to old level, it's the new one now
      CALL read_3D_latbc_from_file(patch_3d, ocean_latbc_timelevel, operators_coeff) ! overwrite old
      CALL save_global(jstep, ocean_latbc_timelevel, patch_2d)
    END DO

    ! calculate the interpolation coefficient
    current_ocean_latbc_dtime_mtime = current_time - ocean_latbc_prev_datetime
    current_ocean_latbc_dtime = REAL(getTotalSecondsTimeDelta(current_ocean_latbc_dtime_mtime, &
      &                                                       ocean_latbc_prev_datetime))
    latbc_f = current_ocean_latbc_dtime / ocean_latbc_dtime

    IF (get_my_global_mpi_id() == 0) write(*,*) "refresh latbc_f = ", latbc_f

    ! interpolate
    l1 = 3 - ocean_latbc_timelevel
    l2 = ocean_latbc_timelevel
    ocean_latbc_data%to(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%to(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%to(l2,:,:,:) ! 1 for temperature
    ocean_latbc_data%so(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%so(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%so(l2,:,:,:) ! 2 for salinity
    ocean_latbc_data%eta_c(3,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%eta_c(l1,:,:) + &
    &    latbc_f *ocean_latbc_data%eta_c(l2,:,:)
    ocean_latbc_data%vn(3,:,:,:) = &
    & (1-latbc_f)*ocean_latbc_data%vn(l1,:,:,:) + &
    &    latbc_f *ocean_latbc_data%vn(l2,:,:,:)

  ! At this point, apply routines are ready to be applied as needed

  END SUBROUTINE update_ocean_latbc

  SUBROUTINE apply_ocean_ssh_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    INTEGER :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: blockNo, jc, jg
    INTEGER :: start_index_c, end_index_c
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_ssh_latbc"

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)

    ! assign to proper cells
    all_cells => patch_2d%cells%all
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c ! do it later using ordering
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 0 .and. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= 4) THEN ! later do it properly
          ocean_state(jg)%p_prog(mt)%eta_c(jc,blockNo) = & !nnew
          & ocean_latbc_data%eta_c(3,jc,blockNo)
        END IF
        ! nudge cells
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 4 .and. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= 8) THEN ! later do it properly
          latbc_nudge = ocean_latbc_data%nudge_cell(jc,blockNo)
          ocean_state(jg)%p_prog(mt)%eta_c(jc,blockNo) = &
          & (1-latbc_nudge)* ocean_state(jg)%p_prog(mt)%eta_c(jc,blockNo) + &
          &    latbc_nudge * ocean_latbc_data%eta_c(3,jc,blockNo) ! 3 for interpolated layer
        END IF
      END DO ! jc
    END DO ! blockNo

  END SUBROUTINE apply_ocean_ssh_latbc

  SUBROUTINE apply_ocean_velocity_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    INTEGER :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges
    INTEGER :: blockNo, je, jk, jg
    INTEGER :: start_index_e, end_index_e
    INTEGER :: start_level, end_level
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_velocity_latbc"

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)
    start_level = 1

    ! assign to proper edges
    all_edges => patch_2d%edges%all
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
      DO je = start_index_e, end_index_e ! do it later using oredering
        end_level = patch_3d%p_patch_1d(1)%dolic_e(je,blockNo)
        IF (patch_2d%edges%refin_ctrl(je,blockNo) > 0 .and. &
          & patch_2d%edges%refin_ctrl(je,blockNo) <= 8) THEN ! later do it properly
          DO jk = start_level, end_level
            ocean_state(jg)%p_prog(mt)%vn(je,jk,blockNo) = & ! 1 for temperature
            & ocean_latbc_data%vn(3,je,jk,blockNo) ! 3 for interpolated layer
          END DO
        END IF
        ! nudge edges
        IF (patch_2d%edges%refin_ctrl(je,blockNo) > 8 .and. &
          & patch_2d%edges%refin_ctrl(je,blockNo) <= 16) THEN ! later do it properly
          latbc_nudge = ocean_latbc_data%nudge_edge(je,blockNo)
          DO jk = start_level, end_level
            ocean_state(jg)%p_prog(mt)%vn(je,jk,blockNo) = &
            & (1-latbc_nudge)* ocean_state(jg)%p_prog(mt)%vn(je,jk,blockNo) + &
            &    latbc_nudge * ocean_latbc_data%vn(3,je,jk,blockNo) ! 3 for interpolated layer
          END DO
        END IF
      END DO ! je
    END DO ! blockNo

  END SUBROUTINE apply_ocean_velocity_latbc

  SUBROUTINE apply_ocean_tracer_latbc(patch_3d, ocean_state, mt)

    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    INTEGER, INTENT(in) :: mt ! model timelevel, either nold(jg) or nnew(jg)

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: blockNo, jc, jk, jg
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_level, end_level
    REAL(wp) :: latbc_nudge
    CHARACTER(*), PARAMETER :: method_name = module_name//"::apply_ocean_tracer_latbc"

    jg = n_dom ! no support for nested grids!
    patch_2d => patch_3d%p_patch_2d(1)
    start_level = 1

    ! assign to proper cells
    all_cells => patch_2d%cells%all
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c ! do it later using ordering
        end_level = patch_3d%p_patch_1d(1)%dolic_c(jc,blockNo)
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 0 .and. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= 4) THEN ! later do it properly
          DO jk = start_level, end_level
            ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,1) = & ! 1 for temperature
            & ocean_latbc_data%to(3,jc,jk,blockNo) ! 3 for interpolated layer
            ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,2) = & ! 2 for salinity
            & ocean_latbc_data%so(3,jc,jk,blockNo) ! 2 for salinity
          END DO
        END IF
        ! nudge cells
        IF (patch_2d%cells%refin_ctrl(jc,blockNo) > 4 .and. &
          & patch_2d%cells%refin_ctrl(jc,blockNo) <= 8) THEN ! later do it properly
          latbc_nudge = ocean_latbc_data%nudge_cell(jc,blockNo)
          DO jk = start_level, end_level
            ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,1) = & ! 1 for temperature
            & (1-latbc_nudge)* ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,1) + &
            &    latbc_nudge * ocean_latbc_data%to(3,jc,jk,blockNo) ! 3 for interpolated layer
            ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,2) = & ! 2 for salinity
            & (1-latbc_nudge)* ocean_state(jg)%p_prog(mt)%tracer(jc,jk,blockNo,2) + &
            &    latbc_nudge * ocean_latbc_data%so(3,jc,jk,blockNo) ! 3 for interpolated layer
          END DO
        END IF
      END DO ! jc
    END DO ! blockNo

  END SUBROUTINE apply_ocean_tracer_latbc

  !-------------------------------------------------------------------------
  SUBROUTINE read_3D_latbc_from_file(patch_3d, timelevel, operators_coeff) !, has_missValue, missValue) ! assuming no missValues for now
!!! ADD CLEANER NULLIFICATION OF LAND CELLS, VERTICAL FILLING AND SMOOTHING MAYBE
    TYPE(t_patch_3d),TARGET, INTENT(inout) :: patch_3d
    INTEGER, INTENT(in) :: timelevel
    TYPE(t_operator_coeff), INTENT(in), TARGET :: operators_coeff

    LOGICAL :: has_missValue(6) = (/ .false., .false., .false., .false., .false., .false. /)
    REAL(wp) :: missValue(6) = (/ -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp, -99999999.0_wp /)
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: varNames(6) = (/ "to ", "so ", "zos", "u  ", "v  ", "vn " /)

    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_stream_id) :: stream_id
    INTEGER :: blockNo, idx, level
    INTEGER :: start_cell_index, end_cell_index, start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_cartesian_coordinates), ALLOCATABLE ::cellVelocity_cc(:,:,:) ! XXX
    CHARACTER(*), PARAMETER :: method_name = module_name//':read_3D_latbc_from_file '
    REAL(wp) :: lat, lon
    !-------------------------------------------------------------------------
    patch_2d => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ocean_salinity(:,:,:) = 0.0_wp
!    has_missValue = .false.
!    missValue     = -99999999.0_wp
!!     IF (initial_salinity_type < 200) RETURN ! not analytic salinity
!
!      CALL message(method_name, ': init from file')
!      CALL init_cell_3D_variable_fromFile(patch_3d, variable=ocean_salinity, name="S", &
!        & has_missValue=has_missValue, missValue=missValue)
!
!    ! copy the initial ocean salinity also to the nudging  salinity
!      IF (no_tracer>1 .AND. type_3dimrelax_salt >0) &
!         ocean_nudge%data_3dimRelax_Salt(:,:,:) = ocean_salinity(:,:,:)
!
!    CALL fillVerticallyMissingValues(patch_3d=patch_3d, ocean_tracer=ocean_salinity,&
!      & has_missValue=has_missValue, missValue=missValue)
!
!    CALL dbg_print('init_ocean_salinity', ocean_salinity(:,:,:), &
!      & module_name,  1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !CALL message (method_name, TRIM(varNames(1))//"...")
    ! read temperature, salinity and sea surface height
    CALL openInputFile(stream_id, TRIM(ocean_latbc_filename), patch_2d)

    CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
      & variable_name=TRIM(varNames(1)), fill_array=ocean_latbc_data%to(timelevel,:,:,:), &
      & has_missValue=has_missValue(1),  missValue=missValue(1))
    CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
      & variable_name=TRIM(varNames(2)), fill_array=ocean_latbc_data%so(timelevel,:,:,:), &
      & has_missValue=has_missValue(2),  missValue=missValue(2))
    CALL read_2D_1time(stream_id=stream_id, location=on_cells, &
      & variable_name=TRIM(varNames(3)), fill_array=ocean_latbc_data%eta_c(timelevel,:,:), &
      & has_missValue=has_missValue(3),  missValue=missValue(3))
    IF (ocean_latbc_from_vn) THEN
      CALL read_3D_1time(stream_id=stream_id, location=on_edges, &
        & variable_name=TRIM(varNames(6)), fill_array=ocean_latbc_data%vn(timelevel,:,:,:), &
        & has_missValue=has_missValue(6),  missValue=missValue(6))
    ELSE
      CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
        & variable_name=TRIM(varNames(4)), fill_array=ocean_latbc_data%u(timelevel,:,:,:), &
        & has_missValue=has_missValue(4),  missValue=missValue(4))
      CALL read_3D_1time(stream_id=stream_id, location=on_cells, &
        & variable_name=TRIM(varNames(5)), fill_array=ocean_latbc_data%v(timelevel,:,:,:), &
        & has_missValue=has_missValue(5),  missValue=missValue(5))
    ENDIF

    CALL closeFile(stream_id)

    ! write(0,*) variable
!     CALL work_mpi_barrier()
!!    DO blockNo = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
!      DO idx = start_cell_index, end_cell_index
!        DO level = patch_3d%p_patch_1d(1)%dolic_c(idx,blockNo) + 1, n_zlev
!!           IF (ocean_latbc_data%to(timelevel,idx,level,blockNo)/= 0.0_wp) THEN
!!             CALL warning(method_name, "non-zero variable on land")
!            ocean_latbc_data%to(timelevel,idx,level,blockNo) = 0.0_wp
!            ocean_latbc_data%so(timelevel,idx,level,blockNo) = 0.0_wp
!            ocean_latbc_data%u(timelevel,idx,level,blockNo) = 0.0_wp
!            ocean_latbc_data%v(timelevel,idx,level,blockNo) = 0.0_wp
!!           ENDIF
!        ENDDO
!      ENDDO
!    ENDDO
!    IF (ocean_latbc_from_vn) THEN
!      DO block = all_edges%start_block, all_edges%end_block
!        CALL get_index_range(all_edges, block, start_edge_index, end_edge_index)
!        DO idx = start_edge_index, end_edge_index
!          DO level = patch_3d%p_patch_1d(1)%dolic_e(idx,block) + 1, n_zlev
!!            IF ( variable(idx,level,block) /=  0.0_wp) THEN
!!              CALL warning(method_name, "non-zero variable on land")
!              ocean_latbc_data%vn(timelevel,idx,level,block) = 0.0_wp
!!            ENDIF
!          ENDDO
!        ENDDO
!      ENDDO
!    ENDIF

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    IF (ocean_latbc_from_vn) THEN
      CALL sync_patch_array(sync_e, patch_2d, ocean_latbc_data%vn(timelevel,:,:,:))
    ELSE
      CALL message(method_name, "gvec2cvec ...")
      ALLOCATE(cellVelocity_cc(nproma,n_zlev,patch_2d%alloc_cell_blocks))

      DO blockNo = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, blockNo, start_cell_index, end_cell_index)
        DO idx = start_cell_index, end_cell_index
          lon = patch_2d%cells%center(idx,blockNo)%lon
          lat = patch_2d%cells%center(idx,blockNo)%lat

          DO level = 1, n_zlev
            CALL gvec2cvec (ocean_latbc_data%u(timelevel, idx, level, blockNo), &
                    & ocean_latbc_data%v(timelevel, idx, level, blockNo), &
                    & lon, lat,                         &
                    & cellVelocity_cc(idx, level, blockNo)%x(1), &
                    & cellVelocity_cc(idx, level, blockNo)%x(2), &
                    & cellVelocity_cc(idx, level, blockNo)%x(3), &
                    & patch_2d%geometry_info)

          ENDDO
        ENDDO
      ENDDO

      CALL message(method_name, "map_cell2edges_3D ...")
      CALL map_cell2edges_3D(patch_3d, cellVelocity_cc, ocean_latbc_data%vn(timelevel,:,:,:), operators_coeff)
      CALL sync_patch_array(sync_e, patch_2d, ocean_latbc_data%vn(timelevel,:,:,:))
      CALL message(method_name, "DEALLOCATE ...")
      DEALLOCATE(cellVelocity_cc)
    ENDIF
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%to(timelevel,:,:,:))
    CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%so(timelevel,:,:,:))
    CALL sync_patch_array(sync_c, patch_2d, ocean_latbc_data%eta_c(timelevel,:,:))

!    WRITE(message_text,*) "miss value (", TRIM(varNames(1)), ") = ", missValue(1), has_missValue(1)
!    CALL message(method_name, message_text)
!    WRITE(message_text,*) "miss value (", TRIM(varNames(2)), ") = ", missValue(2), has_missValue(2)
!    CALL message(method_name, message_text)
!    WRITE(message_text,*) "miss value (", TRIM(varNames(3)), ") = ", missValue(3), has_missValue(3)
!    CALL message(method_name, message_text)
!    WRITE(message_text,*) "miss value (", TRIM(varNames(4)), ") = ", missValue(4), has_missValue(4)
!    CALL message(method_name, message_text)
!    WRITE(message_text,*) "miss value (", TRIM(varNames(5)), ") = ", missValue(5), has_missValue(5)
!    CALL message(method_name, message_text)
!    WRITE(message_text,*) "miss value (", TRIM(varNames(6)), ") = ", missValue(6), has_missValue(6)
!    CALL message(method_name, message_text)

  END SUBROUTINE read_3D_latbc_from_file
  !-------------------------------------------------------------------------

  SUBROUTINE save_const_zone_data(patch_2d)
    TYPE(t_patch), POINTER, INTENT(in) :: patch_2d
    INTEGER :: pe
    INTEGER :: jc, je, blockNo
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_index_e, end_index_e
    CHARACTER(3) :: str2
    TYPE(t_subset_range), POINTER :: all_cells, all_edges

    all_cells => patch_2d%cells%all ! try in_domain vs all and see
    all_edges => patch_2d%edges%all ! try in_domain vs all and see
    pe = get_my_global_mpi_id()
    WRITE(str2,'(I3.3)') pe

    ! Global cells
    OPEN(UNIT=131, FILE="L000/c_lon_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=132, FILE="L000/c_lat_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=133, FILE="L000/c_pe_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=134, FILE="L000/c_refin_ctrl_"//str2, STATUS="replace", ACTION="write")
    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c
        WRITE(131,*) patch_2d%cells%center(jc,blockNo)%lon*rad2deg
        WRITE(132,*) patch_2d%cells%center(jc,blockNo)%lat*rad2deg
        WRITE(133,*) pe
        WRITE(134,*) patch_2d%cells%refin_ctrl(jc,blockNo)
      END DO
    END DO
    CLOSE(131)
    CLOSE(132)
    CLOSE(133)
    CLOSE(134)
    ! Global edges
    OPEN(UNIT=135, FILE="L000/e_lon_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=136, FILE="L000/e_lat_"//str2, STATUS="replace", ACTION="write")
    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
      DO je = start_index_e, end_index_e
        WRITE(135,*) patch_2d%edges%center(je,blockNo)%lon*rad2deg
        WRITE(136,*) patch_2d%edges%center(je,blockNo)%lat*rad2deg
      END DO
    END DO
    CLOSE(135)
    CLOSE(136)
  END SUBROUTINE save_const_zone_data

  SUBROUTINE save_global(timestep, timelevel, patch_2d)
    INTEGER, INTENT(IN) :: timestep
    INTEGER, INTENT(IN) :: timelevel
    TYPE(t_patch), POINTER, INTENT(in) :: patch_2d
    INTEGER :: pe
    INTEGER :: jc, je, blockNo
    INTEGER :: start_index_c, end_index_c
    INTEGER :: start_index_e, end_index_e
    CHARACTER(10) :: str1
    CHARACTER(3) :: str2
    CHARACTER(1) :: str3
    TYPE(t_subset_range), POINTER :: all_cells, all_edges

    all_cells => patch_2d%cells%all ! try in_domain vs all and see
    all_edges => patch_2d%edges%all ! try in_domain vs all and see
    pe = get_my_global_mpi_id()
    WRITE(str1,'(I10.10)') timestep
    WRITE(str2,'(I3.3)') pe
    WRITE(str3,'(I1.1)') timelevel

    ! Global cells
    OPEN(UNIT=131, FILE="L000/c_to_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=132, FILE="L000/c_so_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=133, FILE="L000/c_eta_c_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=134, FILE="L000/c_u_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")
    OPEN(UNIT=135, FILE="L000/c_v_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")

    DO blockNo = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, blockNo, start_index_c, end_index_c)
      DO jc = start_index_c, end_index_c
        WRITE(131,*) ocean_latbc_data%to(timelevel,jc,1,blockNo)
        WRITE(132,*) ocean_latbc_data%so(timelevel,jc,1,blockNo)
        WRITE(133,*) ocean_latbc_data%eta_c(timelevel,jc,blockNo)
        WRITE(134,*) ocean_latbc_data%u(timelevel,jc,1,blockNo)
        WRITE(135,*) ocean_latbc_data%v(timelevel,jc,1,blockNo)
      END DO
    END DO
    CLOSE(131)
    CLOSE(132)
    CLOSE(133)
    CLOSE(134)
    CLOSE(135)

    ! Global edges
    OPEN(UNIT=136, FILE="L000/e_vn_"//str3//"_"//str1//"_"//str2, STATUS="replace", ACTION="write")

    DO blockNo = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_index_e, end_index_e)
      DO je = start_index_e, end_index_e
        WRITE(136,*) ocean_latbc_data%vn(timelevel,je,1,blockNo)
      END DO
    END DO
    CLOSE(136)

  END SUBROUTINE save_global
  !--------------------------------------------------------------------------------------
END MODULE mo_ocean_limarea
