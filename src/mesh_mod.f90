module mesh_mod
use halo_mod, only : halo_t
implicit none

type, public :: mesh_t

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: nx                     !global horizontal lentgh
    integer(kind=4) :: panel_ind
    integer(kind=4) :: halo_width

    real(kind=8), allocatable    :: rhx(:,:), rhy(:,:), rhz(:,:)
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

    this%is = is
    this%ie = ie
    this%js = js
    this%je = je
    this%ks = ks
    this%ke = ke

    this%halo_width = halo_width

end subroutine init_mesh

end module mesh_mod
