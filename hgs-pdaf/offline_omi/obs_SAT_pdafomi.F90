!$Id: obs_B_pdafomi.F90 561 2020-11-21 09:59:46Z lnerger $
!> PDAF-OMI observation module for type B observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = B
!!
!! __Observation type B:__
!! The observation type B in this tutorial are the only the observations at
!! the locations (8,5), (12,15), and (4,30). 
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_SAT_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi_obs_f, &
       ONLY: obs_f   ! Declaration of observation data types
  USE PDAFomi_obs_l, &
       ONLY : obs_l

  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_o_sat

  ! Further variables specific for the saturation observations
  CHARACTER(len=100) :: path_obs_sat  = ''      ! Path to saturation observations
  CHARACTER(len=110) :: file_sat_prefix  = ''   ! file name prefix for saturation observations
  CHARACTER(len=110) :: file_sat_suffix  = '.nc'! file name suffix for saturation observations
  CHARACTER(len=110) :: file_syntobs_sat = 'syntobs.nc' ! File name for synthetic observations

  REAL    :: rms_obs_sat      ! Observation error standard deviation
  REAL    :: bias_obs_sat     ! saturation observation bias

  REAL    :: lradius_sat      ! Localization radius for saturation
  REAL    :: sradius_sat      ! Support radius for localization function

  REAL    :: sat_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)
  LOGICAL :: sat_fixed_rmse   ! Whether to use a fixed RMS error or the error provided with the data

  REAL, ALLOCATABLE :: loc_radius_sat(:) ! Localization radius array for saturation



! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_SAT(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assim_pdaf, &
         ONLY: filtertype, use_global_obs, time, istep, state_type, offset
    USE hgsdat, &
         ONLY: nn, x, y, z

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j, s                      ! Counters
    INTEGER :: n_obs_sat
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    REAL(4), ALLOCATABLE :: sat_x(:),sat_y(:),sat_z(:),sat_value(:,:)  ! observation coordinates and values read from the original
    REAL(4), ALLOCATABLE :: std_value(:,:)                         
    INTEGER :: ios
    CHARACTER(len=100) :: sat_file = ''     ! Complete name of head observation file without path
    INTEGER :: stat(100)                 ! Status for NetCDF functions
    INTEGER :: fileid                    ! ID for NetCDF file
    INTEGER :: id_dim, id_x, id_y, id_z, id_sat, id_std, id_obsid
    INTEGER :: i_obs, steps, offset_step
    INTEGER, ALLOCATABLE :: obs_id(:)
    REAL, ALLOCATABLE :: obs_sat_error_p(:) ! PE-local observed head error


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'HGS-PDAF','Assimilate observations - obs type HEAD'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_o_sat) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2
    thisobs%use_global_obs = use_global_obs


! **********************************
! *** Read PE-local observations ***
! **********************************

! Hear we read observation file in netCDF format. Preprocessing routines are in separate files.

   ! Read observations at the current timestep


    sat_file=TRIM(file_sat_prefix)//TRIM(file_sat_suffix)

    s = 1
    stat(s) = NF_OPEN(TRIM(path_obs_sat)//sat_file, NF_NOWRITE, fileid)

    ! Read no. of observations
    s = s + 1
    stat(s) = NF_INQ_DIMID(fileid, 'n_obs', id_dim)
    s = s + 1
    stat(s) = NF_INQ_DIMLEN(fileid, id_dim, n_obs_sat)
    !dim_obs = dim_obs_p

    ! Read total time steps
    s = s + 1
    stat(s) = NF_INQ_DIMID(fileid, 'time', id_dim)
    s = s + 1
    stat(s) = NF_INQ_DIMLEN(fileid, id_dim, steps)


    !write(*,*) thisobs%dim_obs_f
    WRITE(*,*) 'n_obs_sat = ', n_obs_sat
    WRITE(*,*) 'total no. of steps = ', steps

    ALLOCATE(sat_x(n_obs_sat))
    ALLOCATE(sat_y(n_obs_sat))
    ALLOCATE(sat_z(n_obs_sat))
    ALLOCATE(sat_value(n_obs_sat, steps))
    ALLOCATE(std_value(n_obs_sat, steps))
    ALLOCATE(obs_id(n_obs_sat))

    ! Read coordinates
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'x', id_x)
    s = s + 1
    stat(s) = NF_GET_VAR_REAL(fileid, id_x, sat_x)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'y', id_y)
    s = s + 1
    stat(s) = NF_GET_VAR_REAL(fileid, id_y, sat_y)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'z', id_z)
    s = s + 1
    stat(s) = NF_GET_VAR_REAL(fileid, id_z, sat_z)

    ! Read saturation observations and std (all steps)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'Saturation', id_sat)
    s = s + 1
    stat(s) = NF_GET_VAR_REAL(fileid, id_sat, sat_value)
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'std', id_std)
    s = s + 1
    stat(s) = NF_GET_VAR_REAL(fileid, id_std, std_value)
    

    ! Read observations' node index
    s = s + 1
    stat(s) = NF_INQ_VARID(fileid, 'obs_id', id_obsid)
    s = s + 1
    stat(s) = NF_GET_VAR_INT(fileid, id_obsid, obs_id)


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count PE-local number of observations ***

    offset_step = 0
    dim_obs_p = 0
    DO i = 1, n_obs_sat
       IF (ABS(sat_value(i,istep+offset_step)) < 1.01) dim_obs_p=dim_obs_p+1
    ENDDO

    WRITE(*,*) 'valid no. of obs_sat=',dim_obs_p
    ! *** Initialize index vector of observed surface nodes ***
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required

    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

    ! *** Initialize PE-local vectors of observations and std error ***

    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(obs_sat_error_p(dim_obs_p))

    ! *** Initialize coordinate arrays for PE-local observations
    ALLOCATE(ocoord_p(2, dim_obs_p))

    i_obs=0
    DO i = 1, n_obs_sat
       IF (ABS(sat_value(i,istep+offset_step)) < 1.01 .AND. ABS(sat_value(i,istep+offset_step)) > 0.0) THEN
          i_obs = i_obs + 1
          IF ((state_type == 3) .OR. (state_type ==6)) THEN
            thisobs%id_obs_p(1, i_obs) = obs_id(i) + offset(2)
          ELSE
            thisobs%id_obs_p(1, i_obs) = obs_id(i)
          END IF
          obs_p(i_obs) = REAL(sat_value(i,istep+offset_step), 8)
          obs_sat_error_p(i_obs) = REAL(std_value(i,istep+offset_step), 8)
          ocoord_p(1, i_obs) = REAL(sat_x(i),8)
          ocoord_p(2, i_obs) = REAL(sat_y(i),8)
       END IF
    ENDDO

    ALLOCATE(ivar_obs_p(dim_obs_p))

! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    IF (sat_fixed_rmse) THEN

       ! *** Set constant SST observation error ***
       IF (mype_filter == 0) &
            WRITE (*, '(a, 5x, a, f12.3, a)') 'HGS-PDAF', &
            '--- Use global saturation observation error of ', rms_obs_sat, 'm'

       obs_sat_error_p(:) = rms_obs_sat
    ELSE

       ! *** Use variable error from file
       IF (mype_filter == 0) &
            WRITE (*,'(a,5x,a,i7)') 'HGS-PDAF', &
            '--- Use variable saturation observation error from file'
    END IF

  DO i = 1, dim_obs_p
    ivar_obs_p(i) = 1.0 / obs_sat_error_p(i)**2
  ENDDO


! ****************************
! *** De-bias observations ***
! ****************************

    IF (bias_obs_sat/=0.0 .AND. mype_filter == 0) &
         WRITE (*, '(a, 5x, a, f12.3)') &
         'HGS-PDAF', '--- Use global observation bias of ', bias_obs_sat

    obs_p = obs_p - bias_obs_sat


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, lradius_sat, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=11) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_SAT



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_SAT(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       ! observation operator for observed grid point values
       CALL PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    END IF

  END SUBROUTINE obs_op_SAT



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_SAT(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, &   
         ONLY: locweight

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    write(*,*) 'hello from init_dim_obs_l_SAT'
    !CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
    !     locweight, local_range, srange, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_SAT



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE localize_covar_SAT(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, &   
         ONLY: locweight

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


! *************************************
! *** Apply covariance localization ***
! *************************************

    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, lradius_sat, sradius_sat, &
         coords_p, HP_p, HPH)

  END SUBROUTINE localize_covar_SAT

END MODULE obs_SAT_pdafomi
