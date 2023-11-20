module machine_constants
    implicit none
    integer, parameter :: dr=selected_real_kind(p=13,r=300)
    integer, parameter :: di=selected_int_kind(9)
end module

module hgsdat
   use machine_constants
   use iso_c_binding
   implicit none

   character(100) :: filename
   character(100) :: hgspath
   !character(C_CHAR),bind(C,name="lhgsvar_message") :: message(80)
   character(80) :: message
   character(5)  :: insuffix,outsuffix
   character(40) :: prefix
   character(49) :: conc_name
   logical :: isolf, isconc  
 
  ! file id's
   character(30) :: fid_coord
   character(30) :: fid_coordolf
   character(30) :: fid_elem
   character(30) :: fid_h
   character(30) :: fid_holf
   character(30) :: fid_conc
   character(30) :: fid_he
   character(30) :: fid_rn
   character(30) :: fid_ar
   character(30) :: fid_concolf
   character(30) :: fid_sat
   character(30) :: fid_exchflux
   !character(20) 

   integer :: l_prfx
   integer :: status
   integer :: nvar_head

   logical :: hfile,hfilen
   logical :: done

   integer(di) :: nn   ! Number of nodes in pm domain
   integer(di) :: nn_olf  ! Number of nodes in olf domain 
   integer(di) :: ne  ! Number of elements
   integer(di) :: nx, ny, nz
   integer(c_int),bind(C,name="lhgsvar_nn")     :: lhgs_nn
   integer(c_int),bind(C,name="lhgsvar_ne")     :: lhgs_ne
   integer(c_int),bind(C,name="lhgsvar_nnpe")   :: lhgs_nnpe
   integer(c_int),bind(C,name="lhgsvar_nn_olf") :: lhgs_nn_olf
   !integer(c_int) :: hgs_version
   integer :: hgs_version

   real(dr), allocatable,target :: x(:)
   real(dr), allocatable,target :: y(:)
   real(dr), allocatable,target :: z(:)
   real(dr), allocatable,target :: head(:)
   real(dr), allocatable,target :: headolf(:)
   real(kind=4), allocatable,target :: exchflux(:)
   real(kind=4), allocatable,target :: satur(:)
   real(dr), allocatable,target :: conc_He(:)
   real(dr), allocatable,target :: conc_Rn(:)
   real(dr), allocatable,target :: conc_Ar(:)
   real(dr), allocatable,target :: concolf(:)
   real(dr), allocatable,target :: conc(:)
   real(dr), allocatable,target :: param_k(:)

   type(c_ptr),bind(C,name="lhgsvar_head")      :: headptr
   type(c_ptr),bind(C,name="lhgsvar_headolf")   :: headolfptr
   type(c_ptr),bind(C,name="lhgsvar_exchflux")  :: exchfluxptr
   type(c_ptr),bind(C,name="lhgsvar_satur")     :: saturptr
   type(c_ptr),bind(C,name="lhgsvar_conc")      :: concptr
   type(c_ptr),bind(C,name="lhgsvar_concolf")   :: concolfptr
   type(c_ptr),bind(C,name="lhgsvar_x")         :: xptr
   type(c_ptr),bind(C,name="lhgsvar_y")         :: yptr
   type(c_ptr),bind(C,name="lhgsvar_z")         :: zptr


   integer(c_int),allocatable        :: elemincident(:,:)
   integer(c_int),allocatable,target :: elemincident_vec(:)
   integer(c_int),allocatable,target :: elemprop(:)

   type(c_ptr),bind(C,name="lhgsvar_elemincident") :: elemincidentptr
   type(c_ptr),bind(C,name="lhgsvar_elemprop")     :: elempropptr

end module hgsdat
