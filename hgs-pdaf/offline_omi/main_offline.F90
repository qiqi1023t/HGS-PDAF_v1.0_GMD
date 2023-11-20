!$Id$
!BOP
!
! !ROUTINE: main --- Main program for example of PDAF offline implementation
!
! !INTERFACE:
PROGRAM MAIN_OFFLINE

! !DESCRIPTION:
! This is the main program for an example implementation of
! PDAF with domain-decomposition and offline configuration.
!
! In the offline mode, we assume that the ensemble
! integrations are performed in a separate program (model)
! and the forecasted ensemble can be read from files. After
! initializing the ensemble information by reading model
! outputs, a single analysis step is performed. Subsequently,
! the analysis ensemble can be written to files that can be
! used to initialize another ensemble forecast.
!
! Using PDAF for domain-decomposition, the offline
! mode can be used to perform assimilation with domain-
! decomposed models. If the models write results for each
! sub-domain, these can be read here using the same
! parallelization. Then, the filter analysis can be
! performed utilizing this parallelization. If the files
! contain the full model state, PDAF in offline mode
! can be used either on a single processor, or the
! fields can be distributed in this program to utilize
! the parallelization of the filters.
!
! Parameters can be set in the code, or - preferably -
! by command line arguments that are parsed by the
! module PARSER. The format for this is
! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
! The handles are defined in the code before the calls
! to the routine PARSE.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
             init_parallel, finalize_parallel, MPI_INTEGER, &
             MPI_real8
  USE hgsdat
  USE mod_assim_pdaf, &
       ONLY: istep, state_type
  USE mod_assim_hgs_pdaf, &
       ONLY: damp_k, Sr
  USE PDAF_interfaces_module, &   ! PDAF interface definitions
       ONLY: PDAF_set_ens_pointer

  IMPLICIT NONE
!EOP


! Local variables
  INTEGER :: dim_ens, dim_all                 ! Counter
  REAL, POINTER :: sens_pointer(:,:)
  REAL, ALLOCATABLE :: ens_p(:)
  INTEGER :: i  ! counter
  REAL, ALLOCATABLE :: sens_1d(:)
! **********************
! *** Initialize MPI ***
! **********************
 
  CALL init_parallel() ! initializes MPI

! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: IF (mype_world == 0) THEN

     WRITE (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     WRITE (*, '(9x, a)') 'Data assimilation with PDAF'

     IF (npes_world > 1) THEN
        WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     ELSE
        WRITE (*, '(/21x, a/)') 'Running on 1 PE'
     END IF
     WRITE (*, '(/)')

  END IF initscreen


! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***
!   ! subroutine: init_parallel_pdaf(dim_ens, screen)
!   ! Often dim_ens=0 when calling this routine, because the real ensemble size
!   ! is initialized later in the program
!   ! screen is 1 for ouput on screen

  write(*,*) 'Initialize parallel PDAF'
  CALL init_parallel_pdaf(0, 1)

! ------------------------------- The following part is done in init_pdaf ----------------------------
! *** Initialize HGS model information and readout ***
! *** subroutine is: hgs_init(p,instring,outstring,isolf,isconc,version,states_type,obs_type,logging)

  !write(*,*) 'Initialize HGS'
  !!CALL hgs_init('test',001,000,1,1,4)  
  !call hgs_init('r_01/test','001','000',1,1,2,state_types,1,.true.)
       
       ! CALLING SEQUENCE:
       ! Calls: read_coordinates
       ! Calls: read_elements


! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.

  !CALL initialize()
  ! ------------------------------ END ---------------------------------------------------------------

! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  CALL init_pdaf()   ! even though subroutine is stored in init_pdaf_offline.F90
      
       ! CALLING SEQUENCE for user supplied routines:
       ! Calls: init_ens_offline
           ! *** Here states are collected from the HGS output files
           ! Calls: get_pm_heads
           ! Calls: get_olf_heads
           ! Calls: get_satur
           ! Calls: get_k
       ! Calls: init_pdaf_parse
       ! Calls: init_pdaf_info
       ! Calls: PDAF_init


  ! *** Perform analysis ***

  IF (mype_world == 0) &
       WRITE (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

  ! Here we need to add the time step counter

  CALL assimilation_pdaf_offline()
 
      ! CALLING SEQUENCE (example for filtertype(6) = ESTKF)
      ! Calls main programm, i.e., PDAF_put_state_estkf, -estkf_update, -estkf_analysis, ...
          ! CALLING SEQUENCE of USER ROUTINES only
          ! Calls: collect_state_pdaf --> Not 'used' for non-parallelized ESTKF
          ! Calls: init_dim_obs_pdaf
          ! Calls: obs_op_pdaf
          ! Calls: init_obs_pdaf
          ! Calls: prodRinvA_pdaf
          ! Calls: init_obsvar_pdaf
          ! Calls: prepoststep_ens_offline
              ! *** Here updated states are written to new HGS input files
              ! Calls: set_pm_heads
              ! Calls: set_olf_heads
              ! calls: set_satur
  IF (mype_world ==0) THEN
    CALL PDAF_set_ens_pointer(sens_pointer, status)
    dim_ens = size(sens_pointer, 1)
    dim_all = size(sens_pointer)
    ALLOCATE(sens_1d(dim_all))
    sens_1d = reshape(sens_pointer,[dim_all])
  END IF

  CALL MPI_bcast(dim_ens,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPIerr)

  ALLOCATE(ens_p(dim_ens))

  CALL MPI_scatter(sens_1d,dim_ens,MPI_real8, ens_p, dim_ens, MPI_real8, 0, MPI_COMM_WORLD, MPIerr)

  IF (state_type == 1) THEN ! Only head is updated
     head = ens_p
     CALL set_pm_heads()
  ELSE IF (state_type == 2) THEN  ! head+K
     head = ens_p (1:nn)
     ! Here we can check the updated K and constrain the increment
     WRITE(*,*) 'flag_max_increment=',maxval(abs(ens_p (nn+1:nn+ne)-log10(param_k)))
     !CHECK: loop over all elements
     DO i = 1, ne
       !IF ( abs(ens_p(nn+i)-log10(param_k(i))) > 0.1 ) THEN
       !   IF (ens_p(nn+i)-log10(param_k(i)) > 0 ) THEN
       !       ens_p (nn + i) = log10(param_k (i)) + 0.1
       !   ELSE
       !       ens_p (nn + i) = log10(param_k (i)) - 0.1
       !   END IF
       !   !0.1 * (ens_p (nn + i)-param_k(i))
       !END IF
       ens_p (nn + i) = log10(param_k (i)) + damp_k * (ens_p(nn+i) - log10(param_k (i))) 
     END DO       
     param_k = 10.0 ** (ens_p (nn+1 : nn+ne)) !back transform
     CALL set_pm_heads()
     CALL set_k()
  ELSE IF (state_type == 3) THEN  ! head and saturation 
     head = ens_p (1:nn)
     ! Here we check if the udpated saturation is between 0 and 1
     DO i = 1, nn
       IF (ens_p (nn + i) < Sr) THEN
         ens_p (nn + i) = Sr
       ELSE IF (ens_p( nn + i) > 1.0) THEN
         ens_p (nn + i) = 1.0 
       END IF
     END DO
     satur = ens_p (nn+1 : nn+nn)
     CALL set_pm_heads()
     CALL set_satur()
  ELSE IF (state_type == 4) THEN  ! Only saturation upadted
     ! Here we check if the udpated saturation is between 0 and 1
     DO i = 1, nn
       IF (ens_p (i) < Sr) THEN
         ens_p (i) = Sr
       ELSE IF (ens_p(i) > 1.0) THEN
         ens_p (i) = 1.0
       END IF
     END DO
     satur = ens_p     
     CALL set_satur()
  ELSE IF (state_type == 5) THEN  ! saturation and K
     DO i = 1, nn
       IF (ens_p (i) < Sr) THEN
         ens_p (i) = Sr
       ELSE IF (ens_p(i) > 1.0) THEN
         ens_p (i) = 1.0
       END IF
     END DO
     satur = ens_p (1:nn)
     DO i = 1, ne
       ens_p (nn + i) = log10(param_k (i)) + damp_k * (ens_p(nn+i) - log10(param_k (i)))
     END DO
     param_k = 10.0 ** (ens_p (nn+1 : nn+ne))
     CALL set_satur()
     CALL set_k()
  ELSE IF (state_type == 6) THEN !head, saturation and K
     head = ens_p (1:nn)
     DO i = 1, nn
       IF (ens_p (nn + i) < Sr) THEN
         ens_p (nn + i) = Sr
       ELSE IF (ens_p( nn + i) > 1.0) THEN
         ens_p (nn + i) = 1.0
       END IF
     END DO
     WRITE(*,*) 'flag_max_increment_S=',maxval(abs(ens_p (nn+1:nn+nn)-satur))
     satur = ens_p (nn+1 : nn+nn)
     WRITE(*,*) 'flag_max_increment_K=',maxval(abs(ens_p (nn+nn+1:nn+nn+ne)-log10(param_k)))
     !CHECK: loop over all elements
     DO i = 1, ne
       ens_p (nn + nn+ i) = log10(param_k (i)) + damp_k * (ens_p(nn + nn + i) - log10(param_k (i)))
     END DO
     param_k = 10.0 ** ens_p (nn+nn+1:nn+nn+ne)  
     CALL set_pm_heads()
     CALL set_satur()
     CALL set_k() 
  END IF

  ! Synchronize at barrier for exit
  CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
  WRITE (*,*) 'model PE exited: mype ', mype_world


! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

  ! *** Finalize PDAF - print memory and timing information
  CALL finalize_pdaf(0)

  IF (mype_world == 0) &
       WRITE (*, '(/1x, a)') 'PDAF offline mode: END'


! *** Terminate MPI
  CALL finalize_parallel()
 
  IF (mype_world ==0) DEALLOCATE(sens_1d) 
  DEALLOCATE(ens_p)

END PROGRAM MAIN_OFFLINE
