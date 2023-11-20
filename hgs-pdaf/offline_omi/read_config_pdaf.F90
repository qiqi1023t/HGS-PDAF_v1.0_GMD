!$Id: read_config_pdaf.F90 2282 2020-04-16 09:20:44Z lnerger $
!BOP
!
! !ROUTINE: read_config_pdaf - Read configuration for PDAF
!
! !INTERFACE: 
SUBROUTINE read_config_pdaf()  

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.

! !USES:
  USE mod_parallel_pdaf, &                       ! Variables for ensemble parallelization
       ONLY: mype_model, n_modeltasks, task_id
  USE mod_assim_pdaf, &                          ! General variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, dim_bias, &
       screen, screen, step_null, filtertype, subtype, &
       incremental, type_trans, type_sqrt, &
       type_forget, forget, locweight, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       twin_experiment, dim_obs_max, use_global_obs, &
       state_type, ResultPath, istep
  USE mod_assim_hgs_pdaf, &                      ! General variables for assimilation for hgs
       ONLY: delt_obs_head, delt_obs_head_offset, &
       damp_k, Sr
  USE obs_HEAD_pdafomi, &
       ONLY: assim_o_head, path_obs_head, file_head_prefix, &
       file_head_suffix, rms_obs_head, bias_obs_head, &
       lradius_head, sradius_head, head_exclude_diff, &
       head_fixed_rmse
  USE obs_SAT_pdafomi, &
       ONLY: assim_o_sat, path_obs_sat, file_sat_prefix, &
       file_sat_suffix, rms_obs_sat, bias_obs_sat, &
       lradius_sat, sradius_sat, sat_exclude_diff, &
       sat_fixed_rmse

  USE output_pdaf, &                             ! Variables for file output
       ONLY: write_da, write_ens, str_daspec

  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile ='namelist.pdaf'    ! Name of namelist file
  CHARACTER(len=32)  :: handle                      ! Handle for command line parser
  LOGICAL :: printconfig = .TRUE.                   ! Print information on all configuration parameters


  ! General settings
  NAMELIST /pdaf/ n_modeltasks, dim_ens, dim_lag, dim_bias, filtertype, &
       subtype, incremental, type_forget, forget, &
       type_trans, type_sqrt, step_null, locweight, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       twin_experiment, dim_obs_max, use_global_obs, &
       write_da, write_ens, str_daspec, printconfig, &
       twin_experiment, istep

  ! Settings specific for HGS
  NAMELIST /pdaf_hgs/ screen, assim_o_head, path_obs_head, &
       file_head_prefix, file_head_suffix, state_type, &
       rms_obs_head, bias_obs_head, lradius_head, sradius_head, &
       head_exclude_diff, head_fixed_rmse, ResultPath, &
       assim_o_sat, path_obs_sat, file_sat_prefix, &
       file_sat_suffix, rms_obs_sat, bias_obs_sat, &
       lradius_sat, sradius_sat, sat_exclude_diff, &
       sat_fixed_rmse, damp_k, Sr



! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML=pdaf)
  READ (20,NML=pdaf_hgs)
  CLOSE (20)

! *** Add trailing slash to paths ***
  CALL add_slash(path_obs_head)
  CALL add_slash(path_obs_sat)

  !CALL add_slash(path_obs_he)
  CALL add_slash(path_init)

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'HGS-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'HGS-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.2)') 'HGS-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'HGS-PDAF','varscale    ', varscale
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','loctype     ', loctype
     WRITE (*,'(a,5x,a,es10.2)')'HGS-PDAF','lradius_head ', lradius_head
     WRITE (*,'(a,5x,a,es10.2)') 'HGS-PDAF','sradius_head ', sradius_head
     WRITE (*,'(a,5x,a,es10.2)') 'HGS-PDAF','loc_ratio   ', loc_ratio
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','assim_o_head   ', assim_o_head
     WRITE (*,'(a,5x,a,es10.2)')'HGS-PDAF','rms_obs_head ', rms_obs_head
     WRITE (*,'(a,5x,a,es10.2)')'HGS-PDAF','bias_obs_head  ', bias_obs_head
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','head_fixed_rmse', head_fixed_rmse
     WRITE (*,'(a,5x,a,f11.3)') 'HGS-PDAF','head_exclude_diff', head_exclude_diff
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','use_global_obs', use_global_obs
     WRITE (*,'(a,5x,a,i10)')   'HGS-PDAF','dim_lag     ', dim_lag
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','path_obs_head     ', TRIM(path_obs_head)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_head_prefix  ', TRIM(file_head_prefix)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_head_suffix  ', TRIM(file_head_suffix)
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','assim_o_sat   ', assim_o_sat
     WRITE (*,'(a,5x,a,es10.2)')'HGS-PDAF','rms_obs_sat ', rms_obs_sat
     WRITE (*,'(a,5x,a,es10.2)')'HGS-PDAF','bias_obs_sat  ', bias_obs_sat
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','sat_fixed_rmse', sat_fixed_rmse
     WRITE (*,'(a,5x,a,f11.3)') 'HGS-PDAF','sat_exclude_diff', sat_exclude_diff
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','path_obs_sat     ', TRIM(path_obs_sat)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_sat_prefix  ', TRIM(file_sat_prefix)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_sat_suffix  ', TRIM(file_sat_suffix)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','path_init   ', TRIM(path_init)
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_init   ', TRIM(file_init)
     WRITE (*,'(a,5x,a,f10.2)') 'HGS-PDAF','damp_k      ', damp_k
     WRITE (*,'(a,5x,a,f10.2)') 'HGS-PDAF','Sr      ', Sr

     IF (read_inistate) THEN
        WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','file_inistate ', TRIM(file_inistate)
     ENDIF
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','write_da    ', write_da
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','write_ens   ', write_ens
     WRITE (*,'(a,5x,a,a)')     'HGS-PDAF','str_daspec  ',TRIM(str_daspec)
     WRITE (*,'(a,5x,a,l)')     'HGS-PDAF','twin_experiment', twin_experiment
     WRITE (*,'(a,1x,a)') 'HGS-PDAF','-- End of PDAF configuration overview --'

  END IF showconf

END SUBROUTINE read_config_pdaf
! ==============================================================================
!BOP
!
! !ROUTINE: add_slash --- Add trailing slash to path string
!
! !INTERFACE:
SUBROUTINE add_slash(path)

! !DESCRIPTION:
! This routine ensures that a string defining a path
! has a trailing slash.
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=100) :: path  ! String holding the path
!EOP

! *** Local variables ***
  INTEGER :: strlength

! *** Add trailing slash ***
  strlength = LEN_TRIM(path)

  IF (path(strlength:strlength) /= '/') THEN
     path = TRIM(path) // '/'
  END IF

END SUBROUTINE add_slash
