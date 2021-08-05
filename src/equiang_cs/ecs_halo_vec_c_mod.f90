!Object to interpolate vector values to virtual points beyond face edge
!C-staggered grid.

module ecs_halo_vec_c_mod

use halo_mod,          only : halo_t, halo_vec_t
use exchange_halo_mod, only : exchange_t
use ecs_halo_mod,      only : ecs_tile_halo_t

implicit none

type, extends(halo_vec_t) :: ecs_halo_vec_c_t

    integer(kind=4)                        :: ts, te
    class(exchange_t),         allocatable :: exch_halo
    !type(ecs_tile_halo_vec_t), allocatable :: tile(:)

    contains

    procedure :: get_halo_vector => get_ecs_c_vector_halo

end type

contains

subroutine get_ecs_c_vector_halo(this,u,v,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_vec_c_t),  intent(inout) :: this
    class(grid_field_t),      intent(inout) :: u,v
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do_vec(u,v, domain%parcomm)

    !do t=this%ts,this%te
    !    call this%tile(t)%interpv(u%tile(t),v%tile(t),domain%partition%tile(t),halo_width)
    !end do
end subroutine get_ecs_c_vector_halo

end module ecs_halo_vec_c_mod
