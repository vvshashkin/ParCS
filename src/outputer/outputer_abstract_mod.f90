module outputer_abstract_mod

use grid_function_mod, only : grid_function_t
use mesh_mod,              only : mesh_t

implicit none

integer(kind=4), parameter :: std_out_stream = 888

type, abstract, public :: outputer_t

    integer(kind=4) :: out_stream = std_out_stream

contains

    procedure(write_proc), deferred :: write

end type outputer_t

abstract interface
    subroutine write_proc(this, f, mesh, ts, te, file_name, rec_num)
        import outputer_t, grid_function_t, mesh_t
        class(outputer_t),     intent(inout) :: this
        integer(kind=4),       intent(in)    :: ts, te
        type(grid_function_t), intent(inout) :: f(ts:te)
        type(mesh_t),          intent(in)    :: mesh(ts:te)
        character(*),          intent(in)    :: file_name
        integer(kind=4),       intent(in), &
                               optional      :: rec_num
    end subroutine
end interface

contains

end module outputer_abstract_mod
