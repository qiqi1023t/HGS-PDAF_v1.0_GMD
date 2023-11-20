!$Id: init_pdaf_offline.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the offline mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, &
       writepe, mype_filter, abort_parallel, &
       MPI_COMM_WORLD, MPIerr, MPI_real8
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: dim_state_p, screen, filtertype, subtype, dim_ens, &
       incremental, type_forget, forget, dim_lag, &
       rank_analysis_enkf, locweight, &
       type_trans, type_sqrt, delt_obs, &
       state_type, istep, ens_p_1d, offset
  USE obs_HEAD_pdafomi, &            ! Variables for observation type HEAD
       ONLY: assim_o_head, rms_obs_HEAD, lradius_head, sradius_head
  !USE obs_CONC_pdafomi, &            ! Variables for observation type CONC
  !     ONLY: assim_CONC, rms_obs_CONC, file_obs_CONC
  USE obs_SAT_pdafomi, &            ! Variables for observation type SAT
       ONLY: assim_o_sat, rms_obs_SAT
  !USE obs_TRANS_pdafomi, &            ! Variables for observation type TRANS
  !     ONLY: assim_TRANS, rms_obs_TRANS, file_obs_TRANS
  USE hgsdat!, &
      ! ONLY: prefix,insuffix,outsuffix,isolf,isconc,hgs_version
  USE output_pdaf, &
       ONLY: write_da, write_ens, init_output_pdaf


  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
!EOP

! Local variables
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation
  CHARACTER(len=100) :: nmlfile ='namelist.hgs'    ! Name of namelist file
  ! External subroutines
  REAL, ALLOCATABLE :: state_p_tmp(:)
  EXTERNAL :: init_ens_offline  ! Ensemble initialization

  ! Settings specific for HGS
  NAMELIST /hgs/ prefix,insuffix,outsuffix,isolf,isconc,hgs_version

  
! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF - OFFLINE MODE'
  END IF


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 1    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
  state_type = 1   ! variables included in the state vector
  
  dim_ens = n_modeltasks       ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  subtype = 5       ! (5) Offline mode
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs = 1     ! Number of time steps between analysis/assimilation steps

! *** Which observation type to assimilate
  assim_o_head = .true.
  !assim_CONC = .false.
  assim_o_sat = .false.
  !assim_TRANS = .false.

! *** specifications for observations ***
  rms_obs_HEAD = 0.05    ! Observation error standard deviation for observation A
  !rms_obs_CONC = 0.1    ! Observation error standard deviation for observation A
!  rms_obs_SAT = 0.01    ! Observation error standard deviation for observation B
  !rms_obs_TRANS = 0.5    ! Observation error standard deviation for observation C
  
! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  lradius_head = 0  ! Range in grid points for observation domain in local filters
  sradius_head = lradius_head  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting



! *** Read PDAF configuration from namelist ***

  CALL read_config_pdaf()


  write(*,*) 'Starting initialization of HGS'

  ! *** Read namelist file ***
  WRITE(*,*) 'Read HGS namelist file: ',nmlfile

  OPEN (20,file=nmlfile)
  READ (20,NML=hgs)
  CLOSE (20)

  CALL hgs_init()

! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.
  CALL initialize()


! *** Initial Screen output ***
! *** This is optional      ***

  IF (mype_world == 0) call init_pdaf_info()



! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** Here, the full selection of filters is        ***
! *** implemented. In a real implementation, one    ***
! *** reduce this to selected filters.              ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  ALLOCATE (ens_p_1d(dim_state_p * dim_ens))
  ALLOCATE (state_p_tmp(dim_state_p))

  ! Colletc variables in the state vector
  IF (state_type == 1) THEN
    ! state vector includes pm_head only
    CALL get_pm_heads()
    state_p_tmp = head
 
  ELSE IF (state_type == 2) THEN
    ! state vector includes pm_head and K
    CALL get_pm_heads()
    CALL get_k()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+ne) = log10(param_k)
  ELSE IF (state_type == 3) THEN
    ! state vector includes pm_head and saturation 
    CALL get_pm_heads()
    CALL get_satur()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = satur
  ELSE IF (state_type == 4) THEN
   ! state vector includes saturation     
    CALL get_satur()
    state_p_tmp = satur
  ELSE IF (state_type == 5) THEN
   ! state vector includes saturation and K       
    CALL get_satur()
    CALL get_k()
    state_p_tmp (1:nn) = satur
    state_p_tmp (nn+1:nn+ne) = log10(param_k)
  ELSE IF (state_type == 6) THEN
   ! state vector includes pm_head, saturation and K
    CALL get_pm_heads()
    CALL get_satur()
    CALL get_k()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = satur
    state_p_tmp (nn+nn+1:nn+nn+ne) = log10(param_k)
  ELSE IF (state_type == 7) THEN
   ! state vector includes pm_head and tracer -He
    CALL get_pm_heads()
    CALL get_pm_conc_He()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_he
  ELSE IF (state_type == 8) THEN
   ! state vector includes pm_head, tracer -He and K
    CALL get_pm_heads()
    CALL get_pm_conc_He()
    CALL get_k()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_he
    state_p_tmp (nn+nn+1:nn+nn+ne) = log10(param_k)
  ELSE IF (state_type == 9) THEN
   ! state vector includes pm_head and tracer -Rn
    CALL get_pm_heads()
    CALL get_pm_conc_Rn()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_rn
  ELSE IF (state_type == 10) THEN
   ! state vector includes pm_head, tracer -Rn and K          
    CALL get_pm_heads()
    CALL get_pm_conc_Rn()
    CALL get_k()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_rn
    state_p_tmp (nn+nn+1:nn+nn+ne) = log10(param_k)
  ELSE IF (state_type == 11) THEN
   ! state vector includes pm_head and tracer -Ar
    CALL get_pm_heads()
    CALL get_pm_conc_Ar()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_ar
  ELSE IF (state_type == 12) THEN
   ! state vector includes pm_head, tracer -Ar and K
    CALL get_pm_heads()
    CALL get_pm_conc_Ar()
    CALL get_k()
    state_p_tmp (1:nn) = head
    state_p_tmp (nn+1:nn+nn) = conc_ar
    state_p_tmp (nn+nn+1:nn+nn+ne) = log10(param_k)
          


  END IF

    CALL MPI_gather(state_p_tmp, dim_state_p, MPI_real8, ens_p_1d, dim_state_p, MPI_real8, 0, MPI_COMM_WORLD, MPIerr)
    DEALLOCATE(state_p_tmp)

  whichinit: IF (filtertype == 2) THEN
     ! *** EnKF with Monte Carlo init ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Rank of speudo-inverse in analysis
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = 0           ! Smoother lag (not implemented here)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 6,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  ELSE
     ! *** All other filters                       ***
     ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
     filter_param_i(1) = dim_state_p ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = 0           ! Smoother lag (not implemented here)
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, 0, &
          filter_param_i, 7,&
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_ens_offline, &
          screen, status_pdaf)
  END IF whichinit


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

! ******************************
! *** Initialize file output ***
! ******************************

  writepe = .FALSE.
  IF (filterpe) THEN
     IF (mype_filter==0) writepe = .TRUE.
  ENDIF

  IF (write_da) THEN
     IF (istep ==1) CALL  init_output_pdaf(dim_lag, writepe)   ! Initialize Netcdf output
  END IF

  
 DEALLOCATE(ens_p_1d)

END SUBROUTINE init_pdaf
