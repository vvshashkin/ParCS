module config_mod

implicit none

type, public :: config_t
contains
    procedure, public :: parse
end type config_t

contains

subroutine parse(this, config_filename)

    class(config_t),  intent(inout) :: this
    character(len=*), intent(in)    :: config_filename

end subroutine parse

end module config_mod
