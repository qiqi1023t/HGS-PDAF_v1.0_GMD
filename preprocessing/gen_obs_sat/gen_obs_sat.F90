! !Program: obs_gen_head --- generate heads observations for hgs-pdaf from original HGS grok file
!
! !INTERFACE:

PROGRAM obs_gen_sat

! !This program aims to deal with the original HGS grok and observation files.
! From the .grok file we read the observation locations and from the observatiion 
! file we read the observation values. Afterwards we save these information into 
! a netCDF file, which can be later used for for data assimilation. The same is done for 
! std.
 
  IMPLICIT NONE
#include "netcdf.inc"

  INTEGER :: i, j, s, m     ! Counters
  CHARACTER(len=120) :: obs_location_path, obs_location_file, obs_value_path
  CHARACTER(len=120) :: ncfile_in, obs_value_file
  CHARACTER(len=120) :: outpath, outfile, ncfile_out
  CHARACTER(len=150) :: attstr
  INTEGER :: ncid_in, id_dim, dim_x, dim_y, dim_z
  INTEGER :: ncid_out, dimid_time, dimid_nobs
  INTEGER :: id_x, id_y, id_z, id_time, id_head, id_std, id_obsid
  REAL, ALLOCATABLE :: obs (:,:)
  REAL, ALLOCATABLE :: std (:,:)
  INTEGER :: stat(100), stas
  CHARACTER(len=120) :: fileprefix
  CHARACTER(len=120) :: filesuffix
  CHARACTER(len=120) :: prefix, fid_coord, filename
  INTEGER :: dimids(2)
  INTEGER :: day, n_obs, n_layers 
  REAL :: fillvalue, distance_0, distance
  INTEGER :: nn, lhgs_nn, nx, ny, nz
  REAL, ALLOCATABLE :: x(:), y(:), z(:)
  REAL, ALLOCATABLE :: dummy_z(:), id_top(:)
  REAL, ALLOCATABLE :: time(:,:), tmp_id_series(:,:)
  INTEGER, ALLOCATABLE :: obs_id(:)
  INTEGER, ALLOCATABLE :: tmp_id(:)
  REAL, ALLOCATABLE :: xx(:), yy(:), zz(:)
  INTEGER, ALLOCATABLE :: obs_namelist(:)

! *********************
! *** Configuration ***
! *********************

  ! Path to the observation location file
  obs_location_path = '/p/project/cjicg41/jicg4139/input_HGS-PDAF/hgsmodel/01122022/Sequential/Master/HGS/Grokfiles/'

  ! Name of the observation location file
  ! Note: here the file should be preprocessed with vim editor to remove 
  ! 1) the blank lines using :g/^$/d
  ! 2) text lines

  obs_location_file = 'ObservationPoints_selection.dat'  ! For saturation we only use the x and y. z will be calcualted at 1m depth.

  ! Path to the observation value file
  obs_value_path = '/p/project/cjicg41/jicg4139/input_HGS-PDAF/hgsmodel/01122022/Observations/sat/'

  ! Specify observation file
  fileprefix = 'obs'
  filesuffix = '.dat' 
  
  ! Path to and name of output file holding observations
  outpath = './'
  outfile = 'obs_SAT.nc'

  ! HGS coordiantes file
  prefix = 'Flow'
  fid_coord = 'o.coordinates_pm'

  
! *******************************************
! *** Read HEAD data from observation file ***
! *******************************************

  n_obs = 8 
  day = 95
  fillvalue = 1.0E7
  n_layers = 16

  ALLOCATE(obs_namelist(n_obs))

  ! Here we put the name of the observation files

  obs_namelist = (/ 17, 27, 47, 58, 83, 87, 110, 129 /)


  ALLOCATE(x(n_obs))
  ALLOCATE(y(n_obs))
  ALLOCATE(dummy_z(n_obs))
  ALLOCATE(z(n_obs))
  ALLOCATE(time(n_obs,day))
  ALLOCATE(std(n_obs,day))
  ALLOCATE(obs(n_obs,day))

  ! Read x,y,z
  OPEN(20,file=TRIM(obs_location_path)//TRIM(obs_location_file), status='old')
  WRITE(*,*),'open observation location file:',TRIM(obs_location_path)//TRIM(obs_location_file)
  DO i = 1, n_obs
    READ(20,*) x(i),y(i),dummy_z(i)
  END DO
  CLOSE(20)

  ! Read observation values
  ! Note: run the bash script to remove the first line of these files before opening!
  DO i = 1, n_obs
    WRITE(obs_value_file,'(A,I0,A)') TRIM(fileprefix), obs_namelist(i), TRIM(filesuffix)
    OPEN(21,file=TRIM(obs_value_path)//TRIM(obs_value_file), status='old')
    WRITE(*,*),'open observation value file:',TRIM(obs_value_path)//TRIM(obs_value_file)
    DO j = 1,day
      READ(21,*) time(i,j),obs(i,j),std(i,j)
    END DO
    CLOSE(21)
  END DO
 
  ! Find the observation related node index
  ! Read the coordinates file
  filename= trim(prefix)//trim(fid_coord)

  OPEN(11,file = filename, status = 'old', action = 'read', form = 'unformatted', iostat = stas)

  if(stas /= 0) then
     write(*,*) 'FILE ERROR: '//filename
     stop
  endif
  read(11) nn
  lhgs_nn = nn

  allocate(xx(nn),yy(nn),zz(nn))
  xx = 0
  yy = 0
  zz = 0
  read(11) (xx(i),yy(i),zz(i),i=1,nn)

  read(11) nx
  read(11) ny
  read(11) nz
  
  CLOSE(11)


  ! Find the observations' node index

  ! Step 1: find the closet nodes using x and y. 

  WRITE (*,*) 'obs no.','obs_idex','closest distance','x_obs','y_obs','x_node','y_node'

  ALLOCATE(tmp_id(n_obs))
  DO i =1, n_obs
  distance_0 = 999999999.0
    DO j = 1, nn
      distance = sqrt((xx(j)-x(i))*(xx(j)-x(i))+(yy(j)-y(i))*(yy(j)-y(i)))
      IF (distance < distance_0) THEN
        tmp_id(i) = j
        distance_0 = distance
      END IF
    END DO
  WRITE(*,*) i,tmp_id(i),distance_0,x(i),y(i),xx(tmp_id(i)),yy(tmp_id(i))
  END DO


  ! Step 2: find the column nodes with z coordinatesi

  ALLOCATE(tmp_id_series(n_obs,n_layers))
  DO i = 1, n_obs
    m = 0
    DO j = 1, nn
     IF (xx(j) == xx(tmp_id(i)) .AND. yy(j) == yy(tmp_id(i))) THEN
       m = m + 1
       tmp_id_series(i, m) = j
     END IF
    END DO
  END DO

  ! Step 3: find the top layer node

  ALLOCATE(id_top(n_obs))
  DO i = 1, n_obs
    distance_0 = 0.0
    DO j = 1, n_layers
      IF (zz(tmp_id_series(i,j)) > distance_0) THEN
        id_top (i) = tmp_id_series(i,j)
        distance_0 = zz(tmp_id_series(i,j))
      END IF
    END DO
  END DO


  ! Step 4: get the z-coordinate of the observation node
  DO i = 1, n_obs
    z(i) = zz (id_top(i)) - 1.5
  END DO

  ! Step 5: use the x, y and z coordinates of the observation node to search for the closet model node idex

  ALLOCATE (obs_id(n_obs))
  obs_id = 0

  WRITE (*,*) 'obs no.','obs_idex','closest distance','x_obs','y_obs','z_obs','x_node','y_node','z_node'
  DO i = 1, n_obs
  distance_0 = 999999999.0
    DO j = 1, nn
      distance = sqrt((xx(j)-x(i))*(xx(j)-x(i))+(yy(j)-y(i))*(yy(j)-y(i))+(zz(j)-z(i))*(zz(j)-z(i))) 
      IF (distance < distance_0) THEN
        obs_id(i) = j
        distance_0 = distance
      END IF    
    END DO
    WRITE (*,*),i,obs_id(i),distance_0,x(i),y(i),z(i),xx(obs_id(i)),yy(obs_id(i)),zz(obs_id(i))
  END DO 

  ! *** Initialize output file

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  WRITE (*,*) 'Write HEAD observations to file: ',TRIM(ncfile_out)

  s = 1
  stat(s) = NF_CREATE(ncfile_out, 0, ncid_out)

  attstr  = 'HEAD observations'
  s = s + 1
  stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
       TRIM(attstr))
  s = s + 1
  stat(s) = NF_PUT_ATT_REAL(ncid_out, NF_GLOBAL, '_FillValue', NF_REAL, 1, REAL(fillvalue,4)) 

  ! Define dimensions
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'n_obs', n_obs, dimid_nobs)
  s = s + 1
  stat(s) = NF_DEF_DIM(ncid_out, 'time', day, dimid_time)

  ! Define variables
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'x', NF_REAL, 1, dimid_nobs, id_x)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'y', NF_REAL, 1, dimid_nobs, id_y)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'z', NF_REAL, 1, dimid_nobs, id_z)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'obs_id', NF_INT, 1, dimid_nobs, id_obsid)
  s = s + 1


  dimids(1) = dimid_nobs
  dimids(2) = dimid_time
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'time', NF_REAL, 2, dimids, id_time)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'Saturation', NF_REAL, 2, dimids, id_HEAD)
  s = s + 1
  stat(s) = NF_DEF_VAR(ncid_out, 'std', NF_REAL, 2, dimids, id_std)
  s = s + 1
  stat(s) = NF_ENDDEF(ncid_out)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in init of output file, no.', i
  END DO

  ! Write x,y,z
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_x, REAL(x, 4))
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_y, REAL(y, 4))
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_z, REAL(z, 4))

  ! Write observations' node index
  s = s + 1
  stat(s) = NF_PUT_VAR_INT(ncid_out, id_obsid, obs_id)

  ! Write time
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_time, REAL(time, 4))
  
  ! Write HEAD
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_HEAD, REAL(obs, 4))

  ! Write std
  s = s + 1
  stat(s) = NF_PUT_VAR_REAL(ncid_out, id_std, REAL(std, 4))

  ! Close file
  s = 1
  stat(s) = NF_CLOSE(ncid_out)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in write of output file, no.', i
  END DO

  DEALLOCATE(x,y,z,time,std,obs,dummy_z)

  DEALLOCATE(xx,yy,zz,obs_namelist)
  DEALLOCATE(obs_id, tmp_id ,tmp_id_series, id_top)
  WRITE (*,*) '--- Done ---'

END PROGRAM obs_gen_sat
