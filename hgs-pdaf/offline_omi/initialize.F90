!$Id: initialize.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the 2D offline example for PDAF
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the model information for
! PDAF. Here, the global size of the model domain, the global size
! of the model state vector and the sizes for decomposition of the 
! state vector need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! 2020-11 - OS+WK initialized coupling to HGS
! Define the dimension of state vector
!
! !USES:
  USE mod_assim_pdaf, &        ! model variables for PDAF
        ONLY: dim_state_p, state_type, offset
  USE hgsdat, &                     ! HGS variables
        ONLY: nn, nn_olf, ne
       
  IMPLICIT NONE

! !ARGUMENTS:

!EOP

! *** Model specifications ***
 ! if the model were a regular rectangular mesh, grid and state dim could be initilaized in the following way
 ! nx = XXX    ! Extent of grid in x-direction
 ! ny = YYY    ! Extent of grid in y-direction
 ! nz = ZZZ    ! Extent of grid in z-direction
 ! nn = nx * ny * nz


  IF (state_type == 1) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads only'
    dim_state_p = nn 

  ELSE IF (state_type == 2) THEN
    WRITE (*, '(/14x, a)') 'Update: pm_heads + pm_K'
    dim_state_p = nn + ne    
    
  ELSE IF (state_type == 3) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + saturation'
    dim_state_p = nn + nn 

  ELSE IF (state_type == 4) THEN
    WRITE (*, '(/14x, a)') 'Update: Sturation only'
    dim_state_p = nn 

  ELSE IF (state_type == 5) THEN
    WRITE (*, '(/14x, a)') 'Update: Saturation + K'
    dim_state_p = nn + ne

  ELSE IF (state_type == 6) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + saturation + K'
    dim_state_p = nn + nn + ne
  
  ELSE IF (state_type == 7) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He'
    dim_state_p = nn + nn

  ELSE IF (state_type == 8) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He + K'
    dim_state_p = nn + nn + ne

  ELSE IF (state_type == 9) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Rn'
    dim_state_p = nn + nn

  ELSE IF (state_type == 10) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Rn + K'
    dim_state_p = nn + nn + ne

  ELSE IF (state_type == 11) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Ar'
    dim_state_p = nn + nn

  ELSE IF (state_type == 12) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Ar + K'
    dim_state_p = nn + nn + ne
  
  ELSE IF (state_type == 13) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Ar'
    dim_state_p = nn + nn + nn

  ELSE IF (state_type == 14) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Ar + K'
    dim_state_p = nn + nn + nn + ne

  ELSE IF (state_type == 15) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Rn'
    dim_state_p = nn + nn + nn

  ELSE IF (state_type == 16) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Rn + K'
    dim_state_p = nn + nn + nn + ne

  ELSE IF (state_type == 17) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Ar+Rn'
    dim_state_p = nn + nn + nn

  ELSE IF (state_type == 18) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -Ar+Rn + K'
    dim_state_p = nn + nn + nn + ne

  ELSE IF (state_type == 19) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Rn+Ar'
    dim_state_p = nn + nn + nn

  ELSE IF (state_type == 20) THEN
    WRITE (*, '(/14x, a)') 'Update: pm heads + tracers -He+Rn+Ar + K'
    dim_state_p = nn + nn + nn + ne

  END IF 


! *** Specify offset of fields in state vector ***
  ALLOCATE(offset(4))

  offset(1)  = 0                           
  offset(2)  = nn
  offset(3)  = nn + ne
  offset(4)  = nn + nn

! *** Screen output ***
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(5x, a, i7)') &
          'Global model state dimension:', dim_state_p


END SUBROUTINE initialize
