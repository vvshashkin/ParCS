module mesh_mod
implicit none

type, public :: mesh_t
    type(tile_mesh_t), allocatable :: tile(:)
    real(kind=8) :: scale !projection scale factor
    real(kind=8) :: omega !angular speed
    real(kind=8), allocatable :: rotation_axis(:) !angular velocity direction
    integer(kind=4) :: ts, te
contains
end type mesh_t

type, public :: tile_mesh_t

    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: nx, ny                     !global horizontal grid dimension
!    integer(kind=4) :: panel_ind
    integer(kind=4) :: halo_width

    real(kind=8), allocatable    :: rx(:,:), ry(:,:), rz(:,:) !cartesian coordinates of mesh points
    real(kind=8), allocatable    :: a1(:,:,:), a2(:,:,:)      !cartesian coordinates of covariant vecs at mesh points
    real(kind=8), allocatable    :: b1(:,:,:), b2(:,:,:)      !cartesian coordinates of contravariant vecs at mesh points
    real(kind=8), allocatable    :: Q(:,:,:)                  !metric tensor at mesh-points
    real(kind=8), allocatable    :: QI(:,:,:)                 !inverse metric tensor at mesh-points
    real(kind=8), allocatable    :: G(:,:)                    !sqrt of metric tensor det at mesh-points
    real(kind=8)                 :: hx, hy                    !horizontal grid step
    real(kind=8)                 :: shift_i, shift_j          !determines shift of the first grid point from the boundary
    real(kind=8)                 :: alpha_0, beta_0           !determines coord start

contains

    procedure, public :: init => init_tile_mesh
    procedure, public :: get_alpha, get_beta
end type tile_mesh_t

contains

subroutine init_tile_mesh(this, is, ie, js, je, ks, ke, halo_width)

    class(tile_mesh_t), intent(out) :: this
    integer(kind=4),    intent(in)  :: is, ie, js, je, ks, ke
    integer(kind=4),    intent(in)  :: halo_width

    integer(kind=4) :: isv, iev, jsv, jev

    isv = is - halo_width; iev = ie + halo_width
    jsv = js - halo_width; jev = je + halo_width

    allocate(this%rx(isv:iev,jsv:jev))
    allocate(this%ry(isv:iev,jsv:jev))
    allocate(this%rz(isv:iev,jsv:jev))
    allocate(this%a1(3,isv:iev,jsv:jev))
    allocate(this%a2(3,isv:iev,jsv:jev))
    allocate(this%b1(3,isv:iev,jsv:jev))
    allocate(this%b2(3,isv:iev,jsv:jev))
    allocate(this%Q (3,isv:iev,jsv:jev)) !3 elements of 2x2 matrix are stored due to symmetricity
    allocate(this%QI(3,isv:iev,jsv:jev))! -'-'-
    allocate(this%G (isv:iev,jsv:jev))

    this%is = is
    this%ie = ie
    this%js = js
    this%je = je
    this%ks = ks
    this%ke = ke

    this%halo_width = halo_width

end subroutine init_tile_mesh

pure function get_alpha(this, i) result(alpha)

    class(tile_mesh_t), intent(in) :: this
    integer(kind=4),    intent(in) :: i

    real(kind=8) :: alpha

    alpha = this%alpha_0 + (i-1+this%shift_i)*this%hx

end function get_alpha

pure function get_beta(this, j) result(beta)

    class(tile_mesh_t), intent(in) :: this
    integer(kind=4),    intent(in) :: j

    real(kind=8) :: beta

    beta = this%beta_0 + (j-1+this%shift_j)*this%hy

end function get_beta

end module mesh_mod
