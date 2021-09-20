module tiles_mod

use parcomm_mod, only : parcomm_global

implicit none

type, public :: tiles_t
    integer(kind=4) :: Ni, Nj, Nk, Nt
    integer(kind=4), allocatable, dimension(:) :: is, ie, js, je, ks, ke
contains
    procedure, public :: init
    procedure, public :: get_points_num
    procedure, public :: check
    procedure, public :: print
end type tiles_t

contains

    subroutine init(this, Nt, Nk, Nj, Ni, is, ie, js, je, ks, ke)
        class(tiles_t),  intent(inout) :: this
        integer(kind=4), intent(in)    :: Nt, Ni, Nj, Nk
        integer(kind=4), intent(in)    :: is(Nt), ie(Nt), js(Nt), je(Nt), ks(Nt), ke(Nt)

        this%Nt = Nt
        this%Nk = Nk
        this%Nj = Nj
        this%Ni = Ni

        this%is = is; this%ie = ie
        this%js = js; this%je = je
        this%ks = ks; this%ke = ke

        call this%check()

    end subroutine init

    subroutine check(this)
        class(tiles_t), intent(in) :: this

        integer(kind=4) :: t
        logical         :: is_test_passed

        is_test_passed = .true.

        do t = 1, this%Nt
            if (this%ie(t)<this%is(t)) is_test_passed = .false.
            if (this%je(t)<this%js(t)) is_test_passed = .false.
            if (this%ke(t)<this%ks(t)) is_test_passed = .false.
        end do

        if (.not.is_test_passed) call parcomm_global%abort('Problem in tiles_mod check')

    end subroutine check

    pure function get_points_num(this, t) result(points_num)

        class(tiles_t),  intent(in) :: this
        integer(kind=4), intent(in) :: t

        integer(kind=4) :: points_num

        points_num = (this%ke(t) - this%ks(t) + 1)* &
                     (this%je(t) - this%js(t) + 1)* &
                     (this%ie(t) - this%is(t) + 1)

    end function get_points_num

    subroutine print(this, t)
        class(tiles_t),  intent(in) :: this
        integer(kind=4), intent(in) :: t

        character(len=1000) :: is, ie, js, je
        write(is,*) this%is(t); write(ie,*) this%ie(t)
        write(js,*) this%js(t); write(je,*) this%je(t)
        print*, ''
        print '(2(1x,a, a), /, 2(1x,a, a))', 'is = ', trim(adjustl(is)), 'ie = ', trim(adjustl(ie)), &
                                             'js = ', trim(adjustl(js)), 'je = ', trim(adjustl(je))

    end subroutine print

end module tiles_mod
