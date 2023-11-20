!$Id: mod_assim_pdaf.F90 2249 2020-04-06 06:42:37Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assim_pdaf

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! Implementation for the 2D online example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! *** This module holds the variables that are not specific to observations ***

! ! Settings for time stepping - available as namelist read-in
  INTEGER :: step_null = 0 ! initial time step of assimilation
  INTEGER :: pdafstep = 1

! ! Settings for observations - available as command line options
  LOGICAL :: use_global_obs ! Whether to use global full obs, of full obs limited to process domains
  LOGICAL :: twin_experiment = .false.   ! Whether to perform a twin experiment with synthetic observations
  INTEGER :: dim_obs_max   ! Expect max. number of observations for synthetic obs.

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
  INTEGER :: subtype      ! Subtype of filter algorithm
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

  INTEGER :: delt_obs = 1
  INTEGER :: delt_obs_offset ! time step offset until first analysis step

  INTEGER :: istep        ! real time step for HGS and PDAF

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
  REAL    :: varscale=1.0 ! Scaling factor for initial ensemble variance
!    ! SEIK/ETKF/LSEIK/ETKFS
  INTEGER :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T

!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                           !   (0) symmetric square root, (1) Cholesky decomposition
!    ! Localization - LSEIK/LETKF/LESTKF
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRADIUS
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance

  INTEGER :: rank_analysis_enkf
  INTEGER :: dim_state              ! Global size of model state
  INTEGER :: dim_state_p            ! PE-local size of model state
  INTEGER :: istep_asml             ! Time step at end of an forecast phase
  LOGICAL :: flag_final=.false.     ! Whether the current is the final analysis step

!    ! Specific for local filters
  INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector

!    ! Variables for adaptive localization radius
  REAL, ALLOCATABLE :: eff_dim_obs(:)            ! Effective observation dimension
  INTEGER :: loctype       ! Type of localization
                    !   (0) Fixed radius defined by lradius
                    !   (1) Variable radius for constant effective observation dimension
  REAL :: loc_ratio        ! Choose lradius so the effective observation dim. is loc_ratio times dim_ens
 

!    ! File output and input - available as as namelist read-in
  LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
  CHARACTER(len=100) :: path_init = '.'         ! Path to initialization files
  CHARACTER(len=110) :: file_init = '.'    ! netcdf file holding distributed initial
                                                ! state and covariance matrix (added is _XX.nc)
  CHARACTER(len=110) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                                ! state (added is _XX.nc)
  CHARACTER(len=200) :: ResultPath = '.'

!    ! Other variables - _NOT_ available as command line options!
  REAL    :: time          ! model time
  INTEGER, ALLOCATABLE :: offset(:) ! Offsets of fields in state vector
  REAL :: coords_l(2)      ! Coordinates of local analysis domain

  INTEGER :: state_type          
        ! new variable for definition of state vector dim used by PDAF
        ! '1' for pm_heads only
        ! '2' for pm_heads+satur
        ! '3' for pm_heads+olf_heads
        ! '4' for pm_heads+olf_heads+satur
        ! '5' for pm_heads+conc
        ! '6' for pm_heads+satur+conc
        ! '7' for pm_heads+olf_heads+conc+olf_conc
        ! '8' for pm_heads+olf_heads+satur+conc+olf_conc
        ! '9' for pm_heads+K
        ! '10' for pm_heads+satur+K
        ! '11' for pm_heads+olf_heads+K
        ! '12' for pm_heads+olf_heads+satur+K

  INTEGER :: obs_type
        ! new variable for definition of observations to be used by PDAF
        ! '1' for pm_heads only
        ! '2' for pm_heads+satur
        ! '3' for pm_heads+transpiration
        ! '4' for pm_heads+sat+transpiration
        ! '5' for pm_heads+discharge
        ! '6' for pm_heads+sat+discharge
        ! '7' for pm_heads+sat+transpiration+discharge
        ! '8' for pm_heads+conc
        ! '9' for pm_heads+satur+conc
        ! '10' for pm_heads+transpiration+conc
        ! '11' for pm_heads+sat+transpiration+conc
        ! '12' for pm_heads+discharge+conc
        ! '13' for pm_heads+sat+discharge+conc
        ! '14' for pm_heads+sat+transpiration+discharge+conc

REAL, ALLOCATABLE :: ens_p_1d(:)



END MODULE mod_assim_pdaf
