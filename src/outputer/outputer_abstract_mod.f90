module outputer_abstract_mod

use grid_field_mod, only : grid_field_t
use domain_mod,     only : domain_t

implicit none

integer(kind=4), parameter :: std_out_stream = 888

type, abstract, public :: outputer_t

    integer(kind=4) :: out_stream = std_out_stream

contains

    procedure(write_proc), deferred :: write
    procedure, public  :: write_horizontal_vec => write_horizontal_vec

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

subroutine write_horizontal_vec(this,u,v,domain,file_name,rec_num)

    use parcomm_mod,  only : parcomm_global

    class(outputer_t),   intent(inout) :: this
    type(grid_field_t),  intent(inout) :: u, v
    type(domain_t),      intent(in)    :: domain
    character(*),        intent(in)    :: file_name
    integer(kind=4),     intent(in)    :: rec_num

    call parcomm_global%abort("write_horizontal_vector is not implemented for this outputer")

end subroutine write_horizontal_vec

end module outputer_abstract_mod
