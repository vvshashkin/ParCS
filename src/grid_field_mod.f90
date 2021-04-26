module grid_field_mod

implicit none

type, public :: grid_field_t
    integer(kind=4)            :: bs, be ! block_t-array bounds
    type(block_t), allocatable :: block(:)
contains
    ! procedure, public :: init => grid_field_init
end type grid_field_t

type, public :: block_t
    real(kind=8), allocatable :: p(:,:,:)
    integer(kind=4)            :: panel_ind ! index of grid panel hosted the block
    integer(kind=4)            :: is, ie, js, je, ks, ke, nvi, nvj, nvk ! p-array bounds
contains
    procedure, public :: init => block_init
end type block_t

contains

subroutine block_init(this, panel_ind, is, ie, js, je, ks, ke, nvi, nvj, nvk)

    class(block_t),  intent(out) :: this
    integer(kind=4), intent(in)  :: panel_ind
    integer(kind=4), intent(in)  :: is, ie, js, je, ks, ke, nvi, nvj, nvk

    if (is>ie .or. js>je .or. ks>ke .or. nvi<0 .or. nvj<0 .or. nvk<0) then
        print*, 'Error! Problem with grid function initialization! Abort!'
        stop
    end if

    allocate( this%p(is-nvi : ie+nvi, js-nvj : je+nvj, ks-nvk : ke+nvk) )

    this%panel_ind = panel_ind
    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

    this%nvi = nvi
    this%nvj = nvj
    this%nvk = nvk

end subroutine block_init

end module grid_field_mod
