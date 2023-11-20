!$Id: prepoststep_ens_offline.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_offline --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!! 
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances. 
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE hgsdat, &
    ONLY: nn, nn_olf, prefix, outsuffix, head, satur, headolf, conc, concolf
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, writepe
  USE mod_assim_pdaf, &
       ONLY: step_null,state_type, offset, istep, istep_asml
  USE output_pdaf, &
       ONLY: write_da, write_netcdf_pdaf, write_netcdf_pdaf_ens, &
       write_pos_da, write_ens, write_pos_da_ens


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (not relevant for offline mode)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seik_update    (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Calls: PDAF_add_increment
! Calls: PDAF_seik_TtimesA
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: dgesv (LAPACK)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, member, domain      ! counters
  !INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
  !LOGICAL, SAVE :: firstio = .TRUE.    ! File output is peformed for first time?
  LOGICAL, SAVE :: firsttime = .TRUE.  ! Routine is called for first time?
  REAL :: invdim_ens                   ! Inverse ensemble size
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  REAL, ALLOCATABLE :: rmse(:)                ! estimated RMS error
  REAL, ALLOCATABLE :: var_p(:)     ! model state variances
  REAL, ALLOCATABLE :: field(:,:)      ! global model field
  CHARACTER(len=3) :: ensstr           ! String for ensemble member
  CHARACTER(len=2) :: stepstr         ! String for time step
  CHARACTER(len=3) :: anastr          ! String for call type (initial, forecast, analysis)
  CHARACTER(len=1) :: typestr       ! Character indicating call type
  

 ! Variables for parallelization - global fields
  !INTEGER :: offset   ! Row-offset according to domain decomposition
  !REAL, ALLOCATABLE :: variance(:)     ! local variance
  !REAL, ALLOCATABLE :: ens(:,:)       ! global ensemble
  !REAL, ALLOCATABLE :: state(:)       ! global state vector
  !REAL,ALLOCATABLE :: ens_p_tmp(:,:) ! Temporary ensemble for some PE-domain
  !REAL,ALLOCATABLE :: state_p_tmp(:) ! Temporary state for some PE-domain

! **********************
! *** INITIALIZATION ***
! **********************

  
  IF (mype_filter == 0) THEN
     
     IF (.not.firsttime) THEN
        WRITE (*,'(a, 8x,a)') 'HGS-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE 
        WRITE (*,'(a, 8x,a)') 'HGS-PDAF', 'Analyze forecast state ensemble'
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF


  ! Allocate fields
  ALLOCATE(var_p(dim_p))

  ! Initialize numbers
  invdim_ens = 1.0 / REAL(dim_ens)
  invdim_ensm1 = 1.0 / REAL(dim_ens-1)

! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! *** Compute mean state
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'HGS-PDAF', '--- compute ensemble mean'

  ! local  
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)

  ! *** Compute sampled variances ***
  var_p(:) = 0.0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        var_p(j) = var_p(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  var_p(:) = invdim_ensm1 * var_p(:)


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  ! Later on if more than one variables are included in the state vector
  ! We should calculate this for each variable separately

  ALLOCATE (rmse(1))
  DO i = 1, dim_p
     rmse = rmse + var_p(i)
  ENDDO
  rmse = SQRT(rmse / dim_p)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to sampled variance: ', rmse


! **************************
! *** Write output files ***
! **************************

  write_pos_da = istep_asml
  write_pos_da_ens = istep_asml 
  output: IF (write_da) THEN
     ! *** Write output to NetCDF files

!     IF ((step - step_null) > 0) THEN
      IF (firsttime) THEN
        ! *** write forecasted state fields ***
        CALL write_netcdf_pdaf('f', write_pos_da, step, dim_p, state_p, 1, rmse, writepe)
      ELSE 
         !*** write assimilated state fields ***
        CALL write_netcdf_pdaf('a', write_pos_da, step, dim_p, state_p, 1, rmse, writepe)

     END IF

    IF (write_ens) THEN
      IF (firsttime) THEN    
!        ! *** write forecasted state fields ***
        CALL write_netcdf_pdaf_ens('f', write_pos_da_ens, step, dim_p, ens_p, 1, rmse, writepe, dim_ens)

      ELSE 
!        ! *** write assimilated state fields ***
        CALL write_netcdf_pdaf_ens('a', write_pos_da_ens, step, dim_p, ens_p, 1, rmse, writepe, dim_ens)
     END IF

    END IF
  DEALLOCATE (rmse)
 ENDIF output


! **********************************************
! *** Compute RMS errors for smoothed states ***
! **********************************************

!  IF (dim_lag > 0 .AND. step > 0) THEN
!     CALL compute_rms_smoother_pdaf(step, dim_lag, dim_p, dim_ens, state_p, var_p)
!  END IF


! *********************************************************
! *** Deallocate observation-related arrays             ***
! *** which were allocate in the INIT_DIMOBS_F routines ***
! *********************************************************

  IF (step > 0) THEN
     CALL deallocate_obs_pdafomi()
  END IF

! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(var_p)
  firsttime = .FALSE.


END SUBROUTINE prepoststep_ens_offline
