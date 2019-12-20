module mesh_mod
use halo_mod, only : halo_t
implicit none

type, public :: mesh_t

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: nx                     !global horizontal grid dimension
    integer(kind=4) :: panel_ind
    integer(kind=4) :: halo_width

    real(kind=8), allocatable    :: rhx(:,:), rhy(:,:), rhz(:,:) !cartesian coordinates of mesh points
    real(kind=8), allocatable    :: acov(:,:,:), bcov(:,:,:)       !cartesian coordinates of covariant vecs at mesh points
    real(kind=8), allocatable    :: actv(:,:,:), bctv(:,:,:)       !cartesian coordinates of contravariant vecs at mesh points
    real(kind=8), allocatable    :: Q(:,:,:)                       !metric tensor at mesh-points
    real(kind=8), allocatable    :: QI(:,:,:)                      !inverse metric tensor at mesh-points
    real(kind=8), allocatable    :: G(:,:)                         !metric tensor det at mesh-points
    real(kind=8)                 :: hx !horizontal grid step
    class(halo_t), allocatable   :: halo

contains

    procedure, public :: init => init_mesh

end type mesh_t

contains

subroutine init_mesh(this, is, ie, js, je, ks, ke, halo_width)

    class(mesh_t),   intent(out) :: this
    integer(kind=4), intent(in)  :: is, ie, js, je, ks, ke
    integer(kind=4), intent(in)  :: halo_width

    allocate(this%rhx(is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%rhy(is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%rhz(is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%acov(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%bcov(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%actv(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%bctv(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width))
    allocate(this%Q(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width)) !3 elements of 2x2 matrix are stored due to symmetricity
    allocate(this%QI(3, is-halo_width : ie+halo_width , js-halo_width : je+halo_width))! -'-'-
    allocate(this%G(is-halo_width : ie+halo_width , js-halo_width : je+halo_width)) 

    this%is = is
    this%ie = ie
    this%js = js
    this%je = je
    this%ks = ks
    this%ke = ke

    this%halo_width = halo_width

end subroutine init_mesh

end module mesh_mod
