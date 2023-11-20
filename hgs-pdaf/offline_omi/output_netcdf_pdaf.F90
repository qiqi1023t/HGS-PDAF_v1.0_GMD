!$Id: output_netcdf_pdaf.F90 2271 2020-04-08 13:04:09Z lnerger $
!BOP
!
! !MODULE:
MODULE output_pdaf

! !DESCRIPTION:
! This modules provides routines to initialize
! NetCDF output files for HGS-PDAF and to write
! output into the files.

! !USES:
  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
  CHARACTER(len=100) :: str_daspec='DA'        ! String to identify assimilation experiment
  CHARACTER(len=1)   :: prec_nc = 's'          ! Precision of NetCDF output
  ! (s) single, (d) double
  INTEGER :: write_pos_da = 1                  ! Counter for next time slice to be written
  INTEGER :: write_pos_da_ens
  LOGICAL :: write_da = .true.                 ! Whether to write output file from assimilation
  LOGICAL :: write_ens = .true.                ! Whether to write output file for each individual ensemble member
!EOP

! Private variables
  LOGICAL, PRIVATE :: debugoutput=.FALSE.   ! Write output for debugging
                                            ! (file contains only last writing)
  INTEGER :: nf_prec                        ! Precision of NetCDF output

CONTAINS
!BOP
!
! !ROUTINE: init_output_pdaf - Initialize NetCDf output file
!
! !INTERFACE:
  SUBROUTINE init_output_pdaf(dim_lag, writepe)

! !USES:
    USE mod_assim_pdaf, &
         ONLY: dim_ens

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP


    IF ( writepe .AND. write_da) THEN
       ! Initialize ocean file
       CALL init_ncfile_flow_pdaf(dim_lag, writepe)
       IF (write_ens) THEN
          CALL init_ncfile_flow_pdaf_ens(dim_lag, writepe, dim_ens)
       END IF
    END IF

  END SUBROUTINE init_output_pdaf
!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ncfile_flow_pdaf - Initialize NetCDf output file for flow fields
!
! !INTERFACE:

  SUBROUTINE init_ncfile_flow_pdaf(dim_lag, writepe)

! !USES:
    USE hgsdat, &
        ONLY: nn,lhgs_ne
    USE mod_assim_pdaf, &
        ONLY: ResultPath
           

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
!EOP

! Local variables
    INTEGER :: i                     ! Counter
    CHARACTER(len=100) :: attstr            ! String to write attributes
    CHARACTER(len=200) :: filename   ! Full name of output file
    INTEGER :: s                            ! auxiliary: status counter
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: fileid
    INTEGER :: DimId_nodes, DimId_elements, DimId_iter, dim1
    INTEGER :: VarId_iter, VarId_time, VarId_asmlstep
    INTEGER :: VarId_ha, VarId_hf, VarId_hi
    INTEGER :: VarId_ka, VarId_kf
    INTEGER :: VarId_sa, VarId_sf
    INTEGER :: member
    INTEGER :: dimarray(2)                  ! auxiliary: array dimension

    

    pe0: IF (writepe) THEN

! Print screen information
       IF (prec_nc == 's') THEN
          WRITE (*, '(/a, 1x, a)') 'HGS-PDAF', 'Initialize assimilation NetCDF flow file - single precision'
          nf_prec = NF_FLOAT
       ELSE
          WRITE (*, '(/a, 1x, a)') 'HGS-PDAF', 'Initialize assimilation NetCDF flow file - double precision'
          nf_prec = NF_DOUBLE
       END IF

! ----- open file and write global attributes

       filename=TRIM(ResultPath)//'hgs-flow.'//TRIM(str_daspec)//'.nc'
       s = 1

       stat(s) = NF_CREATE(filename, 0, fileid)
       s = s + 1

       attstr  = 'Flow Assimilation'
       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr))

       s = s + 1


! ----- DEFINE DIMENSIONS ------------------------------

       stat(s) = NF_DEF_DIM(fileid, 'nodes', nn, DimId_nodes)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'elements', lhgs_ne, DimId_elements)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'one', 1, dim1)
       s = s + 1
       stat(s) = NF_DEF_DIM(fileid, 'iteration', NF_UNLIMITED, DimId_iter)
       s = s + 1

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in flow dimension definitions, no.', i
       END DO

! ----- DEFINE VARIABLES ---------------------------------

       s = 1

       !- numbers

       stat(s) = NF_DEF_VAR(fileid, 'iter', NF_INT, 1, DimId_iter, VarId_iter)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'time',   nf_prec, 1, DimId_iter, VarId_time)
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'asmlstep', NF_INT, 1, DimId_iter, VarId_asmlstep)
       s = s + 1

! ----- F I E L D S

       !- scalar variables

       ! Head
       dimarray(1) = DimId_nodes
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'head_a', nf_prec, 2, dimarray(1:2), VarId_ha);
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'head_f', nf_prec, 2, dimarray(1:2), VarId_hf);
       s = s + 1

       dimarray(2) = dim1
       stat(s) = NF_DEF_VAR(fileid, 'head_ini', nf_prec, 2, dimarray(1:2), VarId_hi);
       s = s + 1

       ! K
       dimarray(1) = DimId_elements
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'K_a', nf_prec, 2, dimarray(1:2), VarId_ka);
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'K_f', nf_prec, 2, dimarray(1:2), VarId_kf);
       s = s + 1

       ! Saturation
       dimarray(1) = DimId_nodes
       dimarray(2) = DimId_iter

       stat(s) = NF_DEF_VAR(fileid, 'sat_a', nf_prec, 2, dimarray(1:2), VarId_sa);
       s = s + 1
       stat(s) = NF_DEF_VAR(fileid, 'sat_f', nf_prec, 2, dimarray(1:2), VarId_sf);
       s = s + 1

       DO i = 1,  s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in flow variable definition, no.', i
       END DO

       !- RMS errors

       ! This should be done when more variables are included into the state vector

       
! ----- DEFINE ATTRIBUTES -----------------------------------------

       s = 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',    4, 'time');
       s = s + 1
       attstr  = 'day'
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr)
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     9, 'iteration')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_asmlstep, 'long_name',    17, 'Assimilation step')
       s = s + 1

       ! Fields
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf,      'long_name',   32, 'hydraulic head - forecast')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf,      'units',        5, 'meter')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf,      'field',       19, 'head, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf,      'connections', 20, 'triangles, triangles')
       s = s + 1
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha,      'long_name',   32, 'hydraulic head - analysis')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha,      'units',        5, 'meter')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha,      'field',       19, 'head, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha,      'connections', 20, 'triangles, triangles')
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi,      'long_name',   31, 'hydraulic head - initial')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi,      'units',        5, 'meter')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi,      'field',       19, 'head, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi,      'connections', 20, 'triangles, triangles')
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'long_name',   32, 'K - before analysis')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'units',       15, 'meter per day')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'field',       19, 'K, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'connections', 20, 'triangles, triangles')
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'long_name',   32, 'K - analysis')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'units',       15, 'meter per day')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'field',       19, 'K, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'connections', 20, 'triangles, triangles')

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'long_name',   32, 'Saturation - before analysis')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'units',       15, 'none')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'field',       19, 'Saturation, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'connections', 20, 'triangles, triangles')
       s = s + 1

       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'long_name',   32, 'Saturation - analysis')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'units',       15, 'none')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'field',       19, 'Saturation, scalar, series')
       s = s + 1
       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'connections', 20, 'triangles, triangles')


       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in flow file attribute assignments, no.', i
       END DO

       s = 1

       stat(s) = NF_ENDDEF(fileid)
       s = s + 1
       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing flow NetCDF file'
       END IF
    END IF pe0

  END SUBROUTINE init_ncfile_flow_pdaf

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ncfile_flow_pdaf_ens - Initialize NetCDf output file for flow fields for each ensemble member
!
! !INTERFACE:

  SUBROUTINE init_ncfile_flow_pdaf_ens(dim_lag, writepe, dim_ens)

! !USES:
    USE hgsdat, &
        ONLY: nn,lhgs_ne
    USE mod_assim_pdaf, &
        ONLY: ResultPath


    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim_lag              ! Smoother lag
    LOGICAL, INTENT(in) :: writepe
    INTEGER, INTENT(IN) :: dim_ens
!EOP

! Local variables
    INTEGER :: i                     ! Counter
    CHARACTER(len=100) :: attstr            ! String to write attributes
    CHARACTER(len=200) :: filename   ! Full name of output file
    INTEGER :: s                            ! auxiliary: status counter
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: fileid
    INTEGER :: DimId_nodes, DimId_elements, DimId_iter, dim1
    INTEGER :: VarId_iter, VarId_time, VarId_asmlstep
    INTEGER :: VarId_ha, VarId_hf, VarId_hi
    INTEGER :: VarId_ka, VarId_kf
    INTEGER :: VarId_sa, VarId_sf
    INTEGER :: member
    INTEGER :: dimarray(2)                  ! auxiliary: array dimension
    CHARACTER(len=10) :: member_string


    pe0: IF (writepe) THEN

! Print screen information
       IF (prec_nc == 's') THEN
          WRITE (*, '(/a, 1x, a)') 'HGS-PDAF', 'Initialize assimilation NetCDF flow file for each ensemble member - singl precision'
          nf_prec = NF_FLOAT
       ELSE
          WRITE (*, '(/a, 1x, a)') 'HGS-PDAF', 'Initialize assimilation NetCDF flow file for each ensemble member - double precision'
          nf_prec = NF_DOUBLE
       END IF

! ----- open file and write global attributes
       DO member = 1, dim_ens
          WRITE(member_string,'(i4.4)') member
          filename=TRIM(ResultPath)//'hgs-flow.'//TRIM(str_daspec)//'_'//TRIM(member_string)//'.nc'

          s = 1

          stat(s) = NF_CREATE(filename, 0, fileid)
          s = s + 1

          attstr  = 'Flow Assimilation'
          stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
            TRIM(attstr))

          s = s + 1

! ----- DEFINE DIMENSIONS ------------------------------

          stat(s) = NF_DEF_DIM(fileid, 'nodes', nn, DimId_nodes)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'elements', lhgs_ne, DimId_elements)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'one', 1, dim1)
          s = s + 1
          stat(s) = NF_DEF_DIM(fileid, 'iteration', NF_UNLIMITED, DimId_iter)
          s = s + 1

          DO i = 1,  s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in flow dimension definitions, no.', i
          END DO

! ----- DEFINE VARIABLES ---------------------------------

          s = 1

          !- numbers

          stat(s) = NF_DEF_VAR(fileid, 'iter', NF_INT, 1, DimId_iter, VarId_iter)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'time',   nf_prec, 1, DimId_iter, VarId_time)
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'asmlstep', NF_INT, 1, DimId_iter, VarId_asmlstep)
          s = s + 1

! ----- F I E L D S

          !- scalar variables
          ! Head

          dimarray(1) = DimId_nodes
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'head_a', nf_prec, 2, dimarray(1:2), VarId_ha);
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'head_f', nf_prec, 2, dimarray(1:2), VarId_hf);
          s = s + 1

          dimarray(2) = dim1
          stat(s) = NF_DEF_VAR(fileid, 'head_ini', nf_prec, 2, dimarray(1:2), VarId_hi);
          s = s + 1

          ! K
          dimarray(1) = DimId_elements
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'K_a', nf_prec, 2, dimarray(1:2), VarId_ka);
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'K_f', nf_prec, 2, dimarray(1:2), VarId_kf);
          s = s + 1

          ! Saturation

          dimarray(1) = DimId_nodes
          dimarray(2) = DimId_iter

          stat(s) = NF_DEF_VAR(fileid, 'sat_a', nf_prec, 2, dimarray(1:2), VarId_sa);
          s = s + 1
          stat(s) = NF_DEF_VAR(fileid, 'sat_f', nf_prec, 2, dimarray(1:2), VarId_sf);
          s = s + 1   

          DO i = 1,  s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in flow variable definition, no.', i
          END DO

          !- RMS errors

          ! This should be done when more variables are included into the state
          ! vector


! ----- DEFINE ATTRIBUTES -----------------------------------------

          s = 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time, 'long_name', 4, 'time');
          s = s + 1
          attstr  = 'day'
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time, 'units', LEN_TRIM(attstr), attstr)
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter, 'long_name', 9, 'iteration')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_asmlstep, 'long_name', 17, 'Assimilation step')
          s = s + 1

          ! Fields
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf, 'long_name', 32, 'hydraulic head - forecast')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf, 'units', 5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf, 'field', 19, 'head, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hf, 'connections', 20, 'triangles, triangles')
          s = s + 1
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha, 'long_name', 32, 'hydraulic head - analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha, 'units', 5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha, 'field', 19, 'head, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ha, 'connections', 20, 'triangles, triangles')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi, 'long_name', 31, 'hydraulic head - initial')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi, 'units', 5, 'meter')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi, 'field', 19, 'head, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_hi, 'connections', 20, 'triangles, triangles')
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'long_name',   32, 'K - before analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'units',       15, 'meter per day')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'field',       19, 'K, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_kf,      'connections', 20, 'triangles, triangles')
          s = s + 1
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'long_name',   32, 'K - analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'units',       15, 'meter per day')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'field',       19, 'K, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ka,      'connections', 20, 'triangles, triangles')

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'long_name',   32, 'Sat - before analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'units',       15, 'none')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'field',       19, 'Sat, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sf,      'connections', 20, 'triangles, triangles')
          s = s + 1
          s = s + 1

          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'long_name',   32, 'Sat - analysis')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'units',       15, 'none')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'field',       19, 'Sat, scalar, series')
          s = s + 1
          stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_sa,      'connections', 20, 'triangles, triangles')

          DO i = 1, s - 1
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in flow file attribute assignments, no.', i
          END DO

          s = 1

          stat(s) = NF_ENDDEF(fileid)
          s = s + 1
          stat(1) = NF_CLOSE(fileid)

          IF (stat(1) /= NF_NOERR) THEN
             WRITE(*, *) 'NetCDF error in closing flow NetCDF file'
          END IF
       END DO
    END IF pe0

  END SUBROUTINE init_ncfile_flow_pdaf_ens

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_netcdf_pdaf - Write global fields into NetCDF files
!
! !INTERFACE:
  SUBROUTINE write_netcdf_pdaf(writetype, write_pos_da, iteration, dim, state_l, nfields, rmse, writepe)

! ! USES:
    IMPLICIT NONE

! !ARGUMENTS:
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

    ! Write flow fields
    CALL write_nc_flow_pdaf(writetype, write_pos_da, iteration, dim, state_l, &
         nfields, rmse, writepe)

  END SUBROUTINE write_netcdf_pdaf

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_netcdf_pdaf_ens - Write global fields into NetCDF files
!
! !INTERFACE:
  SUBROUTINE write_netcdf_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

! ! USES:
    IMPLICIT NONE

! !ARGUMENTS:
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da         ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
    INTEGER, INTENT(in) :: dim_ens
!EOP


    ! Write ocean fields
    CALL write_nc_flow_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l,&
         nfields, rmse, writepe, dim_ens)

  END SUBROUTINE write_netcdf_pdaf_ens

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_nc_flow_pdaf - Write flow fields into NetCDF file

!
! !INTERFACE:
  SUBROUTINE write_nc_flow_pdaf(writetype, write_pos_da, iteration, dim, state_l,&
       nfields, rmse, writepe)

! !USES:
    USE mod_assim_pdaf, &
         ONLY: offset, istep_asml, eff_dim_obs, ResultPath, state_type
    USE hgsdat, &
        ONLY: nn,lhgs_ne


    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    REAL, INTENT(in) :: state_l(dim)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

! Local variables
    INTEGER :: FileId                    ! Id of Netcdf file
    CHARACTER(len=200) :: filename   ! Full name of output file
    INTEGER :: s, i
    INTEGER :: VarId_iter, VarId_time, VarId_asmlstep
    INTEGER :: VarId_h, VarId_k, VarId_s 
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2)     ! Position index arrays for writing
    REAL, ALLOCATABLE :: head_value(:), k_value(:), sat_value(:)



    filename=TRIM(ResultPath)//'hgs-flow.'//TRIM(str_daspec)//'.nc'
    ! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write initial flow state to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write flow forecast to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write flow analysis to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       END IF
    END IF

    pe0: IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)
       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs

       s = 1

       stat(s) = NF_INQ_VARID(fileid, "iter", VarId_iter)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "time", VarId_time)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "asmlstep", VarId_asmlstep)
       s = s + 1
       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_ini", VarId_h)
          s = s + 1
       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_a", VarId_h)
          stat(s) = NF_INQ_VARID(fileid, "K_a", VarId_k)
          stat(s) = NF_INQ_VARID(fileid, "sat_a", VarId_s)
          s = s + 1
       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_f", VarId_h)
          stat(s) = NF_INQ_VARID(fileid, "K_f", VarId_k)
          stat(s) = NF_INQ_VARID(fileid, "sat_f", VarId_s)
          s = s + 1

       END IF

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO
! ----- WRITE VARIABLES
       s = 1

       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

       IF (writetype/='s') THEN
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_asmlstep, pos1, nmb, write_pos_da)
          s = s + 1
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, pos1, nmb, istep_asml)
          s = s + 1
          IF (prec_nc == 's') THEN
            !  Write single precision
            stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, pos1, nmb, REAL(istep_asml, 4))
          ELSE
            ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_time, pos1, nmb,REAL(istep_asml, 8))
          END IF
          s = s + 1
       END IF
    END IF pe0

! ----- 2D FIELDS : heads

    !----- hydraulic heads
 
    ALLOCATE (head_value(nn))

    IF (state_type == 1 .OR. state_type == 2 .OR. state_type == 3 .OR. state_type == 6) THEN         
     DO i = 1, nn
       head_value(i) = state_l(i + offset(1))
     END DO

     pe0a: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nn, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_h, pos1vec, nmbvec, &
               REAL(head_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_h, pos1vec, nmbvec, &
               head_value(:))
       END IF

       s = s + 1
     END IF pe0a
    END IF
    DEALLOCATE(head_value)

    ! ----- hydraulic conductivity

    ALLOCATE(k_value(lhgs_ne))

    IF (state_type == 2 .OR. state_type == 5 .OR. state_type == 6) THEN
     IF (state_type == 6) THEN
      DO i = 1, lhgs_ne
         k_value(i) = state_l(i + offset(4))
      END DO
     ELSE
      DO i = 1, lhgs_ne
         k_value(i) = state_l(i + offset(2))
      END DO
     END IF
     IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ lhgs_ne, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_k, pos1vec, nmbvec, &
               REAL(k_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_k, pos1vec, nmbvec, &
               k_value(:))
       END IF

       s = s + 1
     END IF 
    END IF

    DEALLOCATE(k_value)

    ! ----- saturation

    ALLOCATE(sat_value(nn))

    IF (state_type > 2 .AND. state_type < 7) THEN
     IF (state_type == 3 .OR. state_type == 6) THEN
      DO i = 1, nn
         sat_value(i) = state_l(i + offset(2))
      END DO
     ELSE
      DO i = 1, nn
         sat_value(i) = state_l(i + offset(0))
      END DO
     END IF
     IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nn, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_s, pos1vec, nmbvec, &
               REAL(sat_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_s, pos1vec, nmbvec, &
               sat_value(:))
       END IF

       s = s + 1
     END IF
    END IF

    DEALLOCATE(sat_value)


! ----- CLOSE THE FILE

       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF

! ----- Clean up

  END SUBROUTINE write_nc_flow_pdaf

!----------------------------------------------------------------------------
!BOP
!
! !ROUTINE: write_nc_oce_pdaf_ens - Write global ocean fields into NetCDF file
!
! !INTERFACE:
  SUBROUTINE write_nc_flow_pdaf_ens(writetype, write_pos_da, iteration, dim, state_l, &
       nfields, rmse, writepe, dim_ens)

! !USES:
    USE mod_assim_pdaf, &
         ONLY: offset, istep_asml, ResultPath, state_type
    USE hgsdat, &
        ONLY: nn,lhgs_ne

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! !ARGUMENTS:
    CHARACTER(len=1), INTENT(in) :: writetype     ! Write (i) initial, (a) assimilated, (f) forecast fields
    INTEGER, INTENT(in) :: write_pos_da           ! Write position
    INTEGER, INTENT(in) :: iteration              ! Current model time step
    INTEGER, INTENT(in) :: dim                    ! Size of state vector
    INTEGER, INTENT(in) :: dim_ens
    REAL, INTENT(in) :: state_l(dim, dim_ens)              ! State vector
    INTEGER, INTENT(in) :: nfields                ! number of fields in state vector
    REAL, INTENT(in) :: rmse(nfields)             ! Array of RMS errors
    LOGICAL, INTENT(in) :: writepe
!EOP

! Local variables
    INTEGER :: FileId                    ! Id of Netcdf file
    INTEGER :: i, n, member, s              ! Counters
    CHARACTER(len=10) :: member_string
    INTEGER :: pos1, nmb                 ! Position index for writing
    INTEGER :: pos1vec(2), nmbvec(2)     ! Position index arrays for writing
    INTEGER :: pos1vec3(3), nmbvec3(3)   ! Position index arrays for writing
    INTEGER :: VarId_iter, VarId_time    ! Id numbers
    INTEGER :: VarId_asmlstep            ! Number of assimilation step
    INTEGER :: VarId_h, VarId_k, VarId_s ! Ids for fields: head, K, sat
    CHARACTER(len=200) :: filename
    INTEGER :: stat(500)                    ! auxiliary: status array
    REAL, ALLOCATABLE :: head_value(:), k_value(:), sat_value(:)

! Allocate local variables

    DO member = 1, dim_ens
      WRITE(member_string,'(i4.4)') member
      filename=TRIM(ResultPath)//'hgs-flow.'//TRIM(str_daspec)//'_'//TRIM(member_string)//'.nc'

! Print screen information
    IF (writepe) THEN
       IF (writetype == 'i') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write initial flow state to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       ELSE IF (writetype== 'f') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write flow forecast to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       ELSE IF (writetype== 'a') THEN
          WRITE (*, '(a, 8x, a, i9, a, i5)') 'HGS-PDAF', 'Write flow analysis to NetCDF at step ', &
               istep_asml, ' position ', write_pos_da
       END IF
    END IF

! Gather full fields on PE 0 and write into files

    pe0: IF (writepe) THEN

! ----- Open Netcdf File
       stat(1) = NF_OPEN(filename, NF_WRITE, FileId)

       IF (stat(1) /= NF_NOERR) STOP 'nc-file error'

! ----- INQUIRE VARIABLE IDs

       s = 1

       stat(s) = NF_INQ_VARID(fileid, "iter", VarId_iter)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "time",      VarId_time)
       s = s + 1
       stat(s) = NF_INQ_VARID(fileid, "asmlstep", VarId_asmlstep)
       s = s + 1
       IF (writetype=='i') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_ini",       VarId_h)
          s = s + 1
       ELSE IF (writetype=='a') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_a", VarId_h)
          stat(s) = NF_INQ_VARID(fileid, "K_a", VarId_k)
          stat(s) = NF_INQ_VARID(fileid, "sat_a", VarId_s)
          s = s + 1
       ELSE IF (writetype=='f') THEN
          stat(s) = NF_INQ_VARID(fileid, "head_f", VarId_h)
          stat(s) = NF_INQ_VARID(fileid, "K_f", VarId_k)
          stat(s) = NF_INQ_VARID(fileid, "sat_f", VarId_s)
          s = s + 1
       END IF

       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error inquiring variable IDs, no.', i
       END DO

! ----- WRITE VARIABLES
       s = 1

       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1 = write_pos_da
       ELSE
          ! Write keeping only a single time instance
          pos1 = 1
       END IF
       nmb  = 1

       IF (writetype/='s') THEN
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_asmlstep, pos1, nmb, write_pos_da)
          s = s + 1
          stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, pos1, nmb, istep_asml)
          s = s + 1
          IF (prec_nc == 's') THEN
          !  ! Write single precision
             stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, pos1, nmb, REAL(istep_asml, 4))
          ELSE
             ! Write double precision
             stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_time, pos1, nmb, REAL(istep_asml, 8))
          END IF
          s = s + 1
       END IF
    END IF pe0

! ----- 2D FIELDS

    !----- hydraulic heads

    ALLOCATE (head_value(nn))

    IF (state_type == 1 .OR. state_type == 2 .OR. state_type == 3 .OR. state_type == 6) THEN
     DO i = 1, nn
       head_value(i) = state_l(i + offset(1), member)
     END DO
     pe0a: IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nn, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_h, pos1vec, nmbvec, &
               REAL(head_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_h, pos1vec, nmbvec, &
               head_value(:))
       END IF
       s = s + 1
       DO i = 1, s - 1
          IF (stat(i) /= NF_NOERR) &
               WRITE(*, *) 'NetCDF error in writing variables, no.', i
       END DO
     END IF pe0a
    END IF
    DEALLOCATE(head_value)

        ! ----- hydraulic conductivity

    ALLOCATE(k_value(lhgs_ne))

    IF (state_type == 2 .OR. state_type == 5 .OR. state_type == 6) THEN
     IF (state_type == 6) THEN
      DO i = 1, lhgs_ne
         k_value(i) = state_l(i + offset(4),member)
      END DO
     ELSE
      DO i = 1, lhgs_ne
         k_value(i) = state_l(i + offset(2),member)
      END DO
     END IF
     IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ lhgs_ne, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_k, pos1vec, nmbvec, &
               REAL(k_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_k, pos1vec, nmbvec, &
               k_value(:))
       END IF

       s = s + 1
     END IF
    END IF

    DEALLOCATE(k_value)

    ! ----- saturation

    ALLOCATE(sat_value(nn))

    IF (state_type > 2 .AND. state_type < 7) THEN
     IF (state_type == 3 .OR. state_type == 6) THEN
      DO i = 1, nn
         sat_value(i) = state_l(i + offset(2),member)
      END DO
     ELSE
      DO i = 1, nn
         sat_value(i) = state_l(i + offset(0),member)
      END DO
     END IF
     IF (writepe) THEN
       IF (.not. DEBUGOUTPUT) THEN
          ! Normal output
          pos1vec = (/ 1, write_pos_da  /)
       ELSE
          ! Write keeping only a single time instance
          pos1vec = (/ 1, 1  /)
       END IF
       nmbvec  = (/ nn, 1 /)

       IF (prec_nc == 's') THEN
          ! Write single precision
          stat(s) = NF_PUT_VARA_REAL(fileid, VarId_s, pos1vec, nmbvec, &
               REAL(sat_value(:),4))
       ELSE
          ! Write double precision
          stat(s) = NF_PUT_VARA_DOUBLE(fileid, VarId_s, pos1vec, nmbvec, &
               sat_value(:))
       END IF

       s = s + 1
     END IF
    END IF

    DEALLOCATE(sat_value)


! ----- CLOSE THE FILE

       stat(1) = NF_CLOSE(fileid)

       IF (stat(1) /= NF_NOERR) THEN
          WRITE(*, *) 'NetCDF error in closing NetCDF file'
       END IF

! ----- Clean up

   END DO

  END SUBROUTINE write_nc_flow_pdaf_ens

END MODULE output_pdaf
