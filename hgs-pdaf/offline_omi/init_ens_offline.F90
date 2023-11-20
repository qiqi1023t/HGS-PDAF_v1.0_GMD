!$Id: init_ens_offline.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init. It has
! to initialize an ensemble of dim\_ens states.
! For the offline mode, the ensemble will be
! typically read-in from files.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
   USE mod_assim_pdaf, &
        ONLY: state_type, ens_p_1d 
   USE hgsdat
!        ONLY: nn, nn_olf, head, headolf, satur, n_k, param_k
  USE mod_parallel_pdaf, &
       ONLY: mype_world !, npes_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
!       MPIerr, MPIstatus


  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag



! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP



! *** local variables ***
  INTEGER :: member  !i, j, col        ! Counters
  CHARACTER(len=3) :: ensstr           ! String for ensemble member

!   INTEGER, SAVE :: allocflag = 0     ! Flag for memory counting
!   REAL, ALLOCATABLE :: ens(:,:)      ! global ensemble array
!   REAL, ALLOCATABLE :: field(:,:)    ! global model field

!   ! variables and arrays for domain decomposition
!   INTEGER :: offset                  ! Row-offset according to domain decomposition
!   INTEGER :: domain                  ! domain counter
!   REAL,ALLOCATABLE :: ens_p_tmp(:,:) ! Temporary ensemble for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens

  
! ********************************
! *** Read ensemble from files ***
! ********************************

 ! WRITE state vector options
  IF (state_type == 1) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads only'
  ELSEIF (state_type == 2) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads and K'
  ELSEIF (state_type == 3) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads and saturation'
  ELSEIF (state_type == 4) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','Saturation only'
  ELSEIF (state_type == 5) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','Saturation and K'
  ELSEIF (state_type == 6) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_head, saturation and K'
  
  ELSEIF (state_type == 9) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads and olf_heads'
  ELSEIF (state_type == 7) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads, olf_heads, pm_concentration and olf_concentration'
  ELSEIF (state_type == 8) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads, olf_heads, saturation, pm_concentration and olf_concentration'
  ELSEIF (state_type == 10) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads, saturation and K'
  ELSEIF (state_type == 11) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads, olf_heads and K'
  ELSEIF (state_type == 12) THEN
    WRITE (*,'(a,8x,a)') 'state vector includes','pm_heads, olf_heads, saturation and K'
  END IF

  WRITE (*, '(9x, a)') '--- Read initial ensemble from files'
  
    ens_p = reshape (ens_p_1d, (/dim_p, dim_ens/))



  !ELSE IF (state_types == 3) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'

  !    CALL get_pm_heads()

  !    ens_p(1:nn, member) = head(:)

  !    CALL get_olf_heads()

  !    ens_p(nn+1:(nn+nn_olf), member) = headolf(:)

      ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 4) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'
  !    write(*,*) prefix

  !    CALL get_pm_heads()

  !    ens_p(1:nn, member) = head(:)

  !    CALL get_olf_heads()

  !    ens_p(nn+1:(nn+nn_olf), member) = headolf(:)

  !    CALL get_satur()

  !    ens_p((nn+nn_olf+1):(nn+nn_olf+nn), member) = satur(:)

      ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 5) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'

  !    CALL get_pm_heads()
      
  !    ens_p(1:nn, member) = head(:)

  !    CALL get_pm_conc()
      
  !    ens_p((nn+1):(nn+nn), member) = conc(:)

  !    ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 6) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'

  !    CALL get_pm_heads()
      
  !    ens_p(1:nn, member) = head(:)

  !    CALL get_satur()

  !    ens_p((nn+1):(nn+nn), member) = satur(:)
      
  !    CALL get_pm_conc()
      
  !    ens_p((nn+nn+1):(nn+nn+nn), member) = conc(:)

  !    ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 7) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'

  !    CALL get_pm_heads()

  !    ens_p(1:nn, member) = head(:)

  !    CALL get_olf_heads()

  !    ens_p(nn+1:(nn+nn_olf), member) = headolf(:)
      
  !    CALL get_pm_conc()
      
  !    ens_p((nn+nn_olf+1):(nn+nn_olf+nn), member) = conc(:)
      
  !    CALL get_olf_conc()
      
  !    ens_p((nn+nn_olf+nn+1):(nn+nn_olf+nn+nn_olf), member) = concolf(:)

  !    ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 8) THEN
  !  DO member = 1, dim_ens
  !    WRITE (ensstr, '(i0.3)') member

  !    prefix = './r_'//TRIM(ensstr)//'/test'

  !    CALL get_pm_heads()

  !    ens_p(1:nn, member) = head(:)

  !    CALL get_olf_heads()

  !    ens_p(nn+1:(nn+nn_olf), member) = headolf(:)

  !    CALL get_satur()

  !    ens_p((nn+nn_olf+1):(nn+nn_olf+nn), member) = satur(:)
      
  !    CALL get_pm_conc()
      
  !    ens_p((nn+nn_olf+nn+1):(nn+nn_olf+nn+nn), member) = conc(:)
      
  !    CALL get_olf_conc()
      
  !    ens_p((nn+nn_olf+nn+nn+1):(nn+nn_olf+nn+nn+nn_olf), member) = concolf(:)

      ! check if read by writing collected state to text file
  !    OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
  !    WRITE (11, 120) ens_p(:, member)
  !    CLOSE(11)
  !  END DO

  !ELSE IF (state_types == 9) THEN
  ! DO member = 1, dim_ens
  !   WRITE (ensstr, '(i0.3)') member

  !   prefix = './r_'//TRIM(ensstr)//'/test'

  !   CALL get_pm_heads()

  !   ens_p(1:nn, member) = head(:)
 
  !   CALL get_k()  !!!!This is not included in hgs_fun yet!!!

  !   ens_p((nn++1):(nn+n_k), member) = param_k(:, member)
!
!     ! check if read by writing collected state to text file
!     OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
!     WRITE (11, 120) ens_p(:, member)
!     CLOSE(11)
!   END DO
!
!  ELSE IF (state_types == 10) THEN
!   DO member = 1, dim_ens
!     WRITE (ensstr, '(i1)') member
!
!     prefix = './r_'//TRIM(ensstr)//'/test'
!
!     CALL get_pm_heads()
!
!     ens_p(1:nn, member) = head(:)
!
!     CALL get_satur()
!
!     ens_p((nn+1):(nn+nn), member) = satur(:)
!
!     ! CALL get_k()
!
!     ! ens_p((nn+nn_olf+nn+1):(nn+nn_olf+nn+n_k), member) = param_k(:, member)
!
!     ! check if read by writing collected state to text file
!     OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
!     WRITE (11, 120) ens_p(:, member)
!     CLOSE(11)
!   END DO
!
!  ELSE IF (state_types == 11) THEN
!   DO member = 1, dim_ens
!     WRITE (ensstr, '(i1)') member
!
!     prefix = './r_'//TRIM(ensstr)//'/test'
!
!     CALL get_pm_heads()
!
!     ens_p(1:nn, member) = head(:)
!
!     CALL get_olf_heads()
!
!     ens_p(nn+1:(nn+nn_olf), member) = headolf(:)
!
!     ! CALL get_k()
!
!     ! ens_p((nn+nn_olf+nn+1):(nn+nn_olf+nn+n_k), member) = param_k(:, member)
!
!     ! check if read by writing collected state to text file
!     OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
!     WRITE (11, 120) ens_p(:, member)
!     CLOSE(11)
!   END DO
!
!  ELSE IF (state_types == 12) THEN
!   DO member = 1, dim_ens
!     WRITE (ensstr, '(i1)') member
!
!     prefix = './r_'//TRIM(ensstr)//'/test'
!
!     CALL get_pm_heads()
!
!     ens_p(1:nn, member) = head(:)
!
!     CALL get_olf_heads()
!
!     ens_p(nn+1:(nn+nn_olf), member) = headolf(:)
!
!     CALL get_satur()
!
!     ens_p((nn+nn_olf+1):(nn+nn_olf+nn), member) = satur(:)
!
!     ! CALL get_k()
!
!     ! ens_p((nn+nn_olf+nn+1):(nn+nn_olf+nn+n_k), member) = param_k(:, member)
!
!     ! check if read by writing collected state to text file
!     OPEN(11, file = './ens_'//TRIM(ensstr)//'_state.txt', status = 'replace')
!     WRITE (11, 120) ens_p(:, member)
!     CLOSE(11)
!   END DO


! ****************************
! *** Distribute substates ***
! ****************************

  ! This is an example how one could distribute ensemble information over multiple processes


!   mype0c: IF (mype_filter == 0) THEN
!      ! *** Initialize and send sub-state on PE 0 ***
! 
!      ! Initialize sub-ensemble for PE 0
!      DO col = 1, dim_ens
!         DO i=1, dim_p
!            ens_p(i, col) = ens(i, col)
!         END DO
!      END DO
! 
!      ! Define offset in state vectors
!      offset = local_dims(1)
! 
!      DO domain = 2, npes_filter
!         ! Initialize sub-ensemble for other PEs and send sub-arrays
! 
!         ! Allocate temporary buffer array
!         ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))
! 
!         ! Initialize MPI buffer for local ensemble
!         DO col = 1, dim_ens
!            DO i = 1, local_dims(domain)
!               ens_p_tmp(i, col) = ens(i + offset, col)
!            END DO
!         END DO
! 
!         ! Send sub-arrays
!         CALL MPI_send(ens_p_tmp, dim_ens * local_dims(domain), &
!              MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
! 
!         DEALLOCATE(ens_p_tmp)
! 
!         ! Increment offset
!         offset = offset + local_dims(domain)
! 
!      END DO
! 
!   ELSE mype0c
!      ! *** Receive ensemble substates on filter-PEs with rank > 0 ***
! 
!      CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
!           0, 1, COMM_filter, MPIstatus, MPIerr)
!      
!   END IF mype0c
!
! ****************
! *** clean up ***
! ****************

!   IF (mype_filter==0) DEALLOCATE(field, ens)

END SUBROUTINE init_ens_offline
