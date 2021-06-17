module grid_field_mod

use mesh_mod, only : mesh_t, tile_mesh_t

implicit none

type, public :: grid_field_t
    integer(kind=4)                 :: ts, te ! tile_field_t-array bounds
    type(tile_field_t), allocatable :: tile(:)
contains
    procedure, public :: update_s1v1 => update_grid_field_s1v1
    generic :: update => update_s1v1

    procedure, public :: assign_s1v1 => assign_grid_field_s1v1
    procedure, public :: assign_s1   => assign_grid_field_s1
    generic :: assign => assign_s1v1, assign_s1
    procedure, public :: copy => copy_grid_field
    procedure, public :: create_similar => create_similar_grid_field
    procedure, public :: algebraic_norm2 => compute_grid_field_algebraic_norm2
    ! procedure, public :: init => grid_field_init
end type grid_field_t

type, public :: tile_field_t
    real(kind=8), allocatable  :: p(:,:,:)
!    integer(kind=4)            :: panel_ind ! index of grid panel hosted the tile_field
    integer(kind=4)            :: is, ie, js, je, ks, ke, nvi, nvj, nvk ! p-array bounds
contains
    procedure, public :: init => tile_field_init

    procedure, public :: update_s1v1 => tile_field_update_s1v1!v = v + s1*v1
    !update -- generic procedure for v = v + ... operations
    generic :: update => update_s1v1

    procedure, public :: assign_s1   => tile_field_assign_s1  !v = s1
    procedure, public :: assign_s1v1 => tile_field_assign_s1v1!v = s1*v1
    !assign -- generic procedure for v = ... operations
    generic :: assign => assign_s1v1, assign_s1

    procedure, public :: algebraic_norm2 => compute_tile_field_algebraic_norm2
end type tile_field_t

contains

function compute_grid_field_algebraic_norm2(this, mesh, parcomm) result(norm2)
        use parcomm_mod, only : parcomm_t
        use mpi

        class(grid_field_t), intent(in)  :: this
        type(mesh_t),        intent(in)  :: mesh
        type(parcomm_t),     intent(in)  :: parcomm
        real(kind=8)                     :: norm2

        real(kind=8) :: local_norm2
        integer(kind=4) :: t
        integer(kind=4) :: ierr

        local_norm2 = 0.0

        do t = mesh%ts, mesh%te
            local_norm2 = local_norm2 + this%tile(t)%algebraic_norm2(mesh%tile(t))
        end do

        call mpi_allreduce(local_norm2, norm2, 1, mpi_double, mpi_sum, parcomm%comm_w, ierr)
        norm2 = sqrt(norm2)
end function compute_grid_field_algebraic_norm2

function compute_tile_field_algebraic_norm2(this, mesh) result(norm2)
        use parcomm_mod, only : parcomm_t

        class(tile_field_t), intent(in)  :: this
        type(tile_mesh_t),   intent(in)  :: mesh
        real(kind=8)                     :: norm2

        integer(kind=4) :: i, j, k

        norm2 = 0.0_8

        do k = mesh%ks, mesh%ke
            do j = mesh%js, mesh%je
                do i = mesh%is, mesh%ie
                    norm2 = norm2 + this%p(i,j,k)**2
                end do
            end do
        end do
end function compute_tile_field_algebraic_norm2

function create_similar_grid_field(this, mesh) result(grid_field)

    class(grid_field_t), intent(in)  :: this
    type(mesh_t),        intent(in)  :: mesh
    type(grid_field_t)               :: grid_field


    integer(kind=4) :: t

    allocate(grid_field%tile(mesh%ts:mesh%te))

    do t = mesh%ts, mesh%te

        call grid_field%tile(t)%init(this%tile(t)%is, this%tile(t)%ie, &
                                     this%tile(t)%js, this%tile(t)%je, &
                                     this%tile(t)%ks, this%tile(t)%ke, &
                                     this%tile(t)%nvi, this%tile(t)%nvj, this%tile(t)%nvk)
    end do

end function create_similar_grid_field

function copy_grid_field(this, mesh) result(grid_field)

    class(grid_field_t), intent(in)  :: this
    type(mesh_t),        intent(in)  :: mesh
    type(grid_field_t)               :: grid_field

    integer(kind=4) :: t, i, j, k

    allocate(grid_field%tile(mesh%ts:mesh%te))

    do t = mesh%ts, mesh%te

        call grid_field%tile(t)%init(this%tile(t)%is, this%tile(t)%ie, &
                                     this%tile(t)%js, this%tile(t)%je, &
                                     this%tile(t)%ks, this%tile(t)%ke, &
                                     this%tile(t)%nvi, this%tile(t)%nvj, this%tile(t)%nvk)

        do k = this%tile(t)%ks,this%tile(t)%ke
            do j = this%tile(t)%js,this%tile(t)%je
                do i = this%tile(t)%is,this%tile(t)%ie
                    grid_field%tile(t)%p(i,j,k) = this%tile(t)%p(i,j,k)
                end do
            end do
        end do

    end do

end function copy_grid_field

subroutine tile_field_init(this, is, ie, js, je, ks, ke, nvi, nvj, nvk)

    class(tile_field_t),  intent(out) :: this
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

end subroutine tile_field_init

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

subroutine tile_field_update_s1v1(this, a1, scalar1, mesh)

    class(tile_field_t), intent(inout) :: this
    type(tile_field_t),  intent(in)    :: a1
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

end subroutine tile_field_update_s1v1
subroutine tile_field_assign_s1v1(this, a1, scalar1, mesh)

    class(tile_field_t), intent(inout) :: this
    type(tile_field_t),  intent(in)    :: a1
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

end subroutine tile_field_assign_s1v1
subroutine tile_field_assign_s1(this, scalar1, mesh)

    class(tile_field_t), intent(inout) :: this
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

end subroutine tile_field_assign_s1
end module grid_field_mod
