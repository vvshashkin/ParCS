module buffer_mod

implicit none

type, public :: buffer_t

    real(kind=8), allocatable :: p(:)

contains

    procedure, public :: init

end type buffer_t

contains

subroutine init(this, n)

    class(buffer_t) :: this
    integer(kind=4), intent(in) :: n

    allocate(this%p(1:n))

end subroutine init

end module buffer_mod
