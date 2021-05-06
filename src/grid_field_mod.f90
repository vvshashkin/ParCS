module grid_field_mod

use mesh_mod, only : mesh_t, tile_mesh_t

implicit none

type, public :: grid_field_t
    integer(kind=4)                 :: ts, te ! tile_array_t-array bounds
    type(tile_array_t), allocatable :: tile(:)
contains
    procedure, public :: update_s1v1 => update_grid_field_s1v1
    generic :: update => update_s1v1

    procedure, public :: assign_s1v1 => assign_grid_field_s1v1
    procedure, public :: assign_s1   => assign_grid_field_s1
    generic :: assign => assign_s1v1, assign_s1
    ! procedure, public :: init => grid_field_init
end type grid_field_t

type, public :: tile_array_t
    real(kind=8), allocatable  :: p(:,:,:)
!    integer(kind=4)            :: panel_ind ! index of grid panel hosted the tile_array
    integer(kind=4)            :: is, ie, js, je, ks, ke, nvi, nvj, nvk ! p-array bounds
contains
    procedure, public :: init => tile_array_init

    procedure, public :: update_s1v1 => tile_array_update_s1v1!v = v + s1*v1
    !update -- generic procedure for v = v + ... operations
    generic :: update => update_s1v1

    procedure, public :: assign_s1   => tile_array_assign_s1  !v = s1
    procedure, public :: assign_s1v1 => tile_array_assign_s1v1!v = s1*v1
    !assign -- generic procedure for v = ... operations
    generic :: assign => assign_s1v1, assign_s1
end type tile_array_t

contains

subroutine tile_array_init(this, is, ie, js, je, ks, ke, nvi, nvj, nvk)

    class(tile_array_t),  intent(out) :: this
    integer(kind=4),      intent(in)  :: is, ie, js, je, ks, ke, nvi, nvj, nvk

    if (is>ie .or. js>je .or. ks>ke .or. nvi<0 .or. nvj<0 .or. nvk<0) then
        print*, 'Error! Problem with grid function initialization! Abort!'
        stop
    end if

    allocate( this%p(is-nvi : ie+nvi, js-nvj : je+nvj, ks-nvk : ke+nvk) )

    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

    this%nvi = nvi
    this%nvj = nvj
    this%nvk = nvk

end subroutine tile_array_init

subroutine update_grid_field_s1v1(this, f1, scalar1, mesh)

    class(grid_field_t), intent(inout) :: this
    type(grid_field_t),  intent(in)    :: f1
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%update(f1%tile(t), scalar1, mesh%tile(t))
    end do

end subroutine update_grid_field_s1v1

subroutine assign_grid_field_s1(this, scalar1, mesh)

    class(grid_field_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(scalar1, mesh%tile(t))
    end do

end subroutine assign_grid_field_s1

subroutine assign_grid_field_s1v1(this, f1, scalar1, mesh)

    class(grid_field_t), intent(inout) :: this
    type(grid_field_t),  intent(in)    :: f1
    real(kind=8),        intent(in)    :: scalar1
    type(mesh_t),        intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call this%tile(t)%assign(f1%tile(t), scalar1, mesh%tile(t))
    end do

end subroutine assign_grid_field_s1v1

subroutine tile_array_update_s1v1(this, a1, scalar1, mesh)

    class(tile_array_t), intent(inout) :: this
    type(tile_array_t),  intent(in)    :: a1
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = this%p(i,j,k) + scalar1*a1%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_array_update_s1v1
subroutine tile_array_assign_s1v1(this, a1, scalar1, mesh)

    class(tile_array_t), intent(inout) :: this
    type(tile_array_t),  intent(in)    :: a1
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1*a1%p(i,j,k)
            end do
        end do
    end do

end subroutine tile_array_assign_s1v1
subroutine tile_array_assign_s1(this, scalar1, mesh)

    class(tile_array_t), intent(inout) :: this
    real(kind=8),        intent(in)    :: scalar1
    type(tile_mesh_t),   intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                this%p(i,j,k) = scalar1
            end do
        end do
    end do

end subroutine tile_array_assign_s1
end module grid_field_mod
