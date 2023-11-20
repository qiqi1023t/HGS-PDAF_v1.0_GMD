!$Id: mod_assim_pdaf.f90 2179 2020-03-18 18:48:06Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assim_hgs_pdaf

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

! *** This module holds the variables specific for observations in the ocean ***

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs_head        ! time step interval between assimilation steps - Ocean
  INTEGER :: delt_obs_head_offset ! time step offset until first analysis step

!    ! File output and input - available as as namelist read-in

  REAL :: damp_k                  ! damping factor for hydraulic conductivity
  REAL :: Sr                      ! Residual saturation 

END MODULE mod_assim_hgs_pdaf

