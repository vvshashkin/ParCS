module grid_function_mod

    implicit none

type, public :: grid_function_t

    real(kind=8), allocatable :: p(:,:,:)
    integer(kind=4)           :: is, ie, js, je, ks, ke
    integer(kind=4)           :: nvi, nvj, nvk

contains

    procedure, public :: init => grid_function_init

end type grid_function_t

contains


subroutine grid_function_init(this, is, ie, js, je, ks, ke, nvi, nvj, nvk)

    class(grid_function_t), intent(out) :: this
    integer(kind=4),        intent(in)  :: is, ie, js, je, ks, ke, nvi, nvj, nvk

    if (is>ie .or. js>je .or. ks>ke .or. nvi<0 .or. nvj<0 .or. nvk<0) then
        print*, 'Error! Problem with grid function initialization! Abort!'
        stop
    end if

    allocate( this.p(is-nvi : ie+nvi, js-nvj : je+nvj, ks-nvk : ke+nvk) )

    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

    this%nvi = nvi
    this%nvj = nvj
    this%nvk = nvk

end subroutine grid_function_init

end module grid_function_mod
