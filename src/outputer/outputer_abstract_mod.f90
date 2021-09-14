module outputer_abstract_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

integer(kind=4), parameter :: std_out_stream = 888

type, abstract, public :: outputer_t

    integer(kind=4) :: out_stream = std_out_stream

contains

    procedure(write_proc), deferred :: write

end type outputer_t

abstract interface
    subroutine write_proc(this, f, domain, file_name, rec_num)
        import outputer_t, grid_field_t, domain_t
        class(outputer_t),   intent(inout) :: this
        type(grid_field_t),  intent(inout) :: f
        type(domain_t),      intent(in)    :: domain
        character(*),        intent(in)    :: file_name
        integer(kind=4),     intent(in), &
                             optional      :: rec_num
    end subroutine
end interface

contains

end module outputer_abstract_mod
