module tile_mod
implicit none

type, public :: tile_t
    integer(kind=4) :: is, ie, js, je, ks, ke
    integer(kind=4) :: panel_number
contains
    procedure, public :: init
    procedure, public :: check
    procedure, public :: print
    procedure, public :: getind
end type tile_t

private

contains

subroutine init(this, is, ie, js, je, ks, ke, panel_number)
    class(tile_t), intent(inout) :: this
    integer(kind=4), intent(in) :: is, ie, js, je, ks, ke, panel_number

    this%is = is; this%ie = ie
    this%js = js; this%je = je
    this%ks = ks; this%ke = ke

    this%panel_number = panel_number

    call this%check()

end subroutine init

subroutine check(this)
    class(tile_t), intent(in) :: this
    logical :: passed

    passed = .true.
    if (this%ie<this%is) passed = .false.
    if (this%je<this%js) passed = .false.
    if (this%ke<this%ks) passed = .false.
    if (this%panel_number<1 .or. this%panel_number>6) passed = .false.

    if (not(passed)) then
        print*, 'Error in tile_mod!!!'
    end if

end subroutine check

subroutine print(this)
    class(tile_t), intent(in) :: this
    character(len=1000) :: is, ie, js, je, panel_number
    write(is,*) this%is; write(ie,*) this%ie
    write(js,*) this%js; write(je,*) this%je
    write(panel_number,*) this%panel_number
    print*, ''
    print*, 'Panel number = ', trim(adjustl(panel_number))
    print '(2(1x,a, a), /, 2(1x,a, a))', 'is = ', trim(adjustl(is)), 'ie = ', trim(adjustl(ie)), &
                                         'js = ', trim(adjustl(js)), 'je = ', trim(adjustl(je))

end subroutine print

subroutine getind(this, is,ie,js,je,ks,ke,panel_number)
    class(tile_t), intent(in) :: this
    integer(kind=4), intent(out), optional :: is, ie, js, je, ks, ke, panel_number

    if(present(is)) is = this%is
    if(present(ie)) ie = this%ie
    if(present(js)) js = this%js
    if(present(je)) je = this%je
    if(present(ks)) ks = this%ks
    if(present(ke)) ke = this%ke
    if(present(panel_number)) panel_number = this%panel_number

end subroutine getind


end module
