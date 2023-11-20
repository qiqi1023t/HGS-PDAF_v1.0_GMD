!------------------------------------------------------------------------
    subroutine hgs_init()! bind(C,name="lhgs_init")
        use iso_C_binding
        use hgsdat
        implicit none

        !obs_types = obs_type

        if(hgs_version.eq.1)then
          fid_coord    = 'o.coordinates'
          fid_coordolf = 'o.coordinates_overland'
          fid_elem     = 'o.elements'
          fid_h        = 'o.head.'
          fid_holf     = 'o.head_overland.'
          fid_conc     = 'o.conc.'
          fid_concolf  = 'o.conc_overland.'
          fid_sat      = 'o.saturation.'
          fid_exchflux = 'o.exchange_flux.'
        else if(hgs_version.eq.2) then
          fid_coord    = 'o.coordinates_pm'
          fid_coordolf = 'o.coordinates_olf'
          fid_elem     = 'o.elements_pm'
          fid_h        = 'o.head_pm.'
          fid_holf     = 'o.head_olf.'
          fid_conc     = 'o.conc_pm.'
          fid_concolf  = 'o.conc_olf.'
          fid_sat      = 'o.sat_pm.'
          fid_exchflux = 'o.ExchFlux_olf.'
        end if

        hgspath='HGS/'
        !call c_f_pointer(c_loc(p),pchar)
        !s_len = index(pchar,c_null_char)-1

        !call c_f_pointer(c_loc(in),inchar)
        !s_len = index(inchar,c_null_char)-1
        !insuffix = trim(inchar(1:s_len))

        !call c_f_pointer(c_loc(out),outchar)
        !s_len = index(outchar,c_null_char)-1
        !outsuffix = trim(outchar(1:s_len))

        call read_coordinates

        call read_elements

        if(isolf) then
            call get_nn_olf
            allocate(headolf(nn_olf))
            headolfptr = c_loc(headolf)
            allocate(exchflux(nn_olf))
            exchfluxptr = c_loc(exchflux)
        endif

        if(isconc) then
            call get_nn_olf
            allocate(concolf(nn_olf))
            concolfptr = c_loc(concolf)
            allocate(conc(nn))
            concptr = c_loc(conc)
        endif

        allocate(satur(nn))
        saturptr = c_loc(satur)
  
        3000 FORMAT (2x,a,T20,a)
        3001 FORMAT (2x,a,T25,i12)
          write(*,*) '------------------------------------------------------------------------'
          write(*,*) '-- hgs_init '
          write(*,*) '------------------------------------------------------------------------'
          write(*,3000) 'file prefix:',prefix
          write(*,3000) 'insuffix:',insuffix
          write(*,3000) 'outsuffix:',outsuffix
          write(*,3001) 'hgs version:',hgs_version
          write(*,3001) 'number of nodes:',nn
          write(*,3001) 'number of surface nodes:',nn_olf
          write(*,3001) 'number of elements:',lhgs_ne
          write(*,*) '------------------------------------------------------------------------'


    end subroutine hgs_init
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine read_coordinates()! bind(C,name="lhgs_read_coordinates")
        !use iso_C_binding
        use hgsdat
       implicit none

        integer :: i

        !filename= trim(prefix)//'o.coordinates'
        filename= trim(hgspath)//trim(prefix)//trim(fid_coord)
       OPEN(11,file = filename, status = 'old', action = 'read', form = 'unformatted', iostat = status)
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//filename
          stop
       endif
       read(11) nn
       lhgs_nn = nn

       allocate(x(nn),y(nn),z(nn),stat=status)
       if(status /= 0) then
          write(*,*) 'ALLOCATION ERROR: Arrays x,y,z'
          stop
       endif
       x = 0
       y = 0
       z = 0
       read(11) (x(i),y(i),z(i),i=1,nn)

       if(hgs_version.eq.1) then
         read(11) nx,ny,nz
       else if(hgs_version.eq.2) then
         read(11) nx
         read(11) ny
         read(11) nz
       end if

       CLOSE(11)

       allocate(head(nn))
       head = 0

       xptr    = c_loc(x)
       yptr    = c_loc(y)
       zptr    = c_loc(z)
       headptr = c_loc(head)

    end subroutine read_coordinates
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine read_elements()
        use hgsdat
       implicit none

       integer :: i,j
       integer :: iter

       !filename= trim(prefix)//'o.elements'
       filename= trim(hgspath)//trim(prefix)//trim(fid_elem)
       OPEN(11,file = filename, status = 'old', action = 'read', form = 'unformatted', iostat = status)
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//filename
          stop
       endif

       if(hgs_version.eq.1) then
         read(11) lhgs_nnpe,lhgs_ne
       else if(hgs_version.eq.2) then
         read(11) lhgs_nnpe
         read(11) lhgs_ne
       end if
       ne=lhgs_ne

       allocate(elemincident(lhgs_nnpe,lhgs_ne),stat=status)
       if(status /= 0) then
          write(*,*) 'ALLOCATION ERROR: elemincident'
          stop
       endif
       allocate(elemincident_vec(lhgs_nnpe*lhgs_ne),stat=status)
       if(status /= 0) then
          write(*,*) 'ALLOCATION ERROR: elemincident'
          stop
       endif

       elemincident=0
       read(11) ((elemincident(j,i),j=1,lhgs_nnpe),i=1,lhgs_ne)

       allocate(elemprop(lhgs_ne),stat=status)
       if(status /= 0) then
          write(*,*) 'ALLOCATION ERROR: elemprop'
          stop
       endif
        
       elemprop=0
       read(11) (elemprop(i),i=1,lhgs_ne)

       close(11)

       iter=0
       do i=1,lhgs_ne
          do j=1,lhgs_nnpe
              iter=iter+1
              elemincident_vec(iter) = elemincident(j,i)
          end do
       end do

       elemincidentptr = c_loc(elemincident_vec)
       elempropptr     = c_loc(elemprop)

       allocate(param_k(ne))

    end subroutine read_elements
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_pm_heads() !bind(C,name="lhgs_get_heads")
        use iso_C_binding
        use hgsdat
       implicit none

       !filename=trim(prefix)//'o.head.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_h)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn,filename,head)
       else
          write(*,*) 'ERROR in reading file ',filename
          stop
       endif

    end subroutine get_pm_heads
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_pm_heads() !bind(C,name="lhgs_set_heads")
       use iso_C_binding
       use hgsdat
       implicit none

       filename=trim(hgspath)//trim(prefix)//trim(fid_h)//outsuffix
       call write_var8(nn,filename,head)
    end subroutine set_pm_heads
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine test_interface() !bind(C,name="lhgs_test_interface")
       use iso_C_binding
       use hgsdat
       implicit none
       integer :: i

       96 FORMAT(F10.6,/)
       write(*,*) message
       write(*,96) (head(i),i=1,nn)
    end subroutine test_interface
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine read_var4(np,fname,var)
       use hgsdat
       implicit none

       integer :: j, np
       character(*) :: fname
       real(kind=4) :: var(np)

       open(12,file=fname,status='old',action = 'read',form='unformatted')
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//fname
          stop
       endif
       read(12) message
       read(12) (var(j),j=1,np)
       close(12)
    end subroutine read_var4
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine write_var4(np,fname,var)
       use hgsdat
       implicit none

       integer :: j, np
       character(*) :: fname
       real(kind=4) :: var(np)

       open(12,file=fname,status='replace',action = 'write',form='unformatted')
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//fname
          stop
       endif
       write(12) message
       write(12) (var(j),j=1,np)
       close(12)
    end subroutine write_var4
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine read_var8(np,fname,var)
       use hgsdat
       implicit none

       integer :: j, np
       character(*) :: fname
       real(dr) :: var(np)

       open(12,file=fname,status='old',action = 'read',form='unformatted')
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//fname
          stop
       endif
       read(12) message
       read(12) (var(j),j=1,np)
       close(12)
    end subroutine read_var8
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine write_var8(np,fname,var)
       use hgsdat
       implicit none

       integer :: j, np
       character(*) :: fname
       real(dr) :: var(np)

       open(12,file=fname,status='replace',action = 'write',form='unformatted')
       if(status /= 0) then
          write(*,*) 'FILE ERROR: '//fname
          stop
       endif
       write(12) message
       write(12) (var(j),j=1,np)
       close(12)
    end subroutine write_var8
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_olf_heads() !bind(C,name="lhgs_get_olfheads")
       use hgsdat
       use iso_C_binding
       implicit none

       !filename=trim(prefix)//'o.head_overland.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_holf)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn_olf,filename,headolf)
       endif

    end subroutine get_olf_heads
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_nn_olf
       use hgsdat
       implicit none

       !OPEN(11,file = trim(prefix)//'o.coordinates_overland', &
       OPEN(11,file = trim(hgspath)//trim(prefix)//trim(fid_coordolf), &
          status = 'old', &
          action = 'read', &
          form = 'unformatted', &
          iostat = status)
       if(status /= 0) then
          !write(*,*) 'FILE ERROR: '//trim(prefix)//'o.coordinates_overland'
          write(*,*) 'FILE ERROR: '//trim(prefix)//trim(fid_coordolf)
          stop
       endif
       read(11) nn_olf
       lhgs_nn_olf = nn_olf

       close(11)
    end subroutine get_nn_olf
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_olf_heads() !bind(C,name="lhgs_set_olfheads")
       use iso_C_binding
       use hgsdat
       implicit none

       !filename=trim(prefix)//'o.head_overland.'//outsuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_holf)//outsuffix
       call write_var8(nn_olf,filename,headolf)

    end subroutine set_olf_heads
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_exchange_flux() !bind(C,name="lhgs_get_exchangeflux")
       use hgsdat
       use iso_C_binding
       implicit none

       !filename=trim(prefix)//'o.exchange_flux.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_exchflux)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var4(nn_olf,filename,exchflux)
       endif

    end subroutine get_exchange_flux
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_satur() !bind(C,name="lhgs_get_satur")
       use hgsdat
       use iso_C_binding
       implicit none

       !filename=trim(prefix)//'o.saturation.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_sat)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var4(nn,filename,satur)
       endif

    end subroutine get_satur
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_satur() !bind(C,name="lhgs_set_satur")
       use hgsdat
       use iso_C_binding
       implicit none

       !filename=trim(prefix)//'o.saturation.'//outsuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_sat)//outsuffix
       call write_var4(nn,filename,satur)

    end subroutine set_satur
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_pm_conc_He() !bind(C,name="lhgs_get_pm_conc")
        use iso_C_binding
        use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc_pm.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_he)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn,filename,conc_he)
       else
          write(*,*) 'ERROR in reading file ',filename
          stop
       endif

    end subroutine get_pm_conc_He
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_pm_conc_Rn() !bind(C,name="lhgs_get_pm_conc")
        use iso_C_binding
        use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc_pm.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_Rn)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn,filename,conc_rn)
       else
          write(*,*) 'ERROR in reading file ',filename
          stop
       endif

    end subroutine get_pm_conc_Rn
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_pm_conc_Ar() !bind(C,name="lhgs_get_pm_conc")
        use iso_C_binding
        use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc_pm.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_ar)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn,filename,conc_ar)
       else
          write(*,*) 'ERROR in reading file ',filename
          stop
       endif

    end subroutine get_pm_conc_Ar
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_olf_conc() !bind(C,name="lhgs_get_pm_conc")
        use iso_C_binding
        use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc_olf.'//insuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_concolf)//insuffix
       inquire(file=filename,exist=hfilen)
       if(hfilen) then
          done=.false.
          call read_var8(nn_olf,filename,concolf)
       else
          write(*,*) 'ERROR in reading file ',filename
          stop
       endif

    end subroutine get_olf_conc
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_pm_conc() !bind(C,name="lhgs_set_conc")
       use iso_C_binding
       use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc.'//outsuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_conc)//trim(conc_name)//'.'//outsuffix
       call write_var8(nn,filename,conc)

    end subroutine set_pm_conc
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_olf_conc() !bind(C,name="lhgs_set_olfconc")
       use iso_C_binding
       use hgsdat
       implicit none

       !filename=trim(prefix)//'o.conc.'//outsuffix
       filename=trim(hgspath)//trim(prefix)//trim(fid_concolf)//outsuffix
       call write_var8(nn_olf,filename,concolf)

    end subroutine set_olf_conc
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine get_k()
       use iso_C_binding
       use hgsdat
       implicit none

       integer :: i
       filename='K.dat'
       open(16,file=TRIM(filename),status='old')
       do i =1, ne
         read(16,*) param_k(i) 
       end do
       close(16)
    end subroutine get_k   
!------------------------------------------------------------------------
!------------------------------------------------------------------------
    subroutine set_k()
       use iso_C_binding
       use hgsdat
       implicit none

       integer :: i
       filename='K.dat'
       open(16,file=TRIM(filename),status='replace',action = 'write')
       do i = 1, ne
         write(16,*) param_k(i)
       end do
       close(16)
    end subroutine set_k
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    subroutine assign_fnvars(p,plen,in,out) bind(C,name="lhgs_assign_hgsdat")
!        use iso_C_binding
!        use hgsdat
!        implicit none
!
!        byte :: p(41)
!        integer :: plen
!        byte :: in(4)
!        byte :: out(4)
!        character*40 :: pchar
!        character*3  :: inchar
!        character*3  :: outchar
!        integer :: i
!
!        write(pchar,'(40a)') (p(i),i=1,plen)
!        write(inchar,'(3a)') (in(i),i=1,3)
!        write(outchar,'(3a)') (out(i),i=1,3)
!
!        prefix    = trim(pchar)
!        insuffix  = trim(inchar)
!        outsuffix = trim(outchar)
!
!    end subroutine assign_fnvars
!------------------------------------------------------------------------


    subroutine read_obsfile_hgs_heads(fn,var,res)
        implicit none
        character(*),intent(in)  :: fn
        integer,intent(in)  :: var
        real,intent(out)  :: res
        integer :: iostatus
        real  :: t,h,p,s,q,h0,p0,q0,x,y,z
        integer :: node


        open(20, file= fn, status='old', action='read',iostat=iostatus)
        if(iostatus /= 0) then
           write(*,*) 'FILE ERROR: '//fn
           stop
        endif

        read(20,*,iostat=iostatus)
        read(20,*,iostat=iostatus)
        read(20,*,iostat=iostatus)
        do
          read(20,*,iostat=iostatus) t,h,p,s,q,h0,p0,q0,x,y,z,node
          if(iostatus/=0) exit
        end do
        write(*,*) t,h,p,s,q,h0,p0,q0,x,y,z,node
        if(var==1) then
          res = h
        elseif (var==2) then
          res = s
        endif

        write(*,*) res

    end subroutine
