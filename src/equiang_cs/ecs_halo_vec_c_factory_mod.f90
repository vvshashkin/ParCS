!Initialization routine for halo_vec object on eq.cubsph grid
!C-type staggering:
module ecs_halo_vec_c_factory_mod
use ecs_halo_mod,       only : ecs_halo_t
use ecs_halo_vec_c_mod, only : ecs_halo_vec_c_t
use parcomm_mod,        only : parcomm_global

implicit none

!integer, parameter :: corner_halo_width = 5!minimum halo-width to compute 2x2 corner-halo-areas

private
public   :: create_ecs_C_vec_halo_procedure

contains

subroutine create_ecs_C_vec_halo_procedure(halo_out,domain,halo_width)
    use halo_mod,               only : halo_vec_t
    use domain_mod,             only : domain_t
    use exchange_factory_mod,   only : create_symmetric_halo_vec_exchange_C

    class(halo_vec_t), allocatable, intent(out) :: halo_out
    class(domain_t),                intent(in)  :: domain
    integer(kind=4),                intent(in)  :: halo_width

    !locals
    type(ecs_halo_vec_c_t), allocatable :: halo
    integer(kind=4)      :: ex_halo_width = 8
    integer(kind=4)      :: ts, te, is,ie, js, je, nh, t
    real(kind=8)         :: hx

    allocate(halo)
    ts = domain%partition%ts
    te = domain%partition%te
    nh = domain%partition%nh
    halo%ts = ts
    halo%te = te
    !allocate(halo%tile(ts:te))

    halo%exch_halo = create_symmetric_halo_vec_exchange_C(domain%partition, &
                              domain%parcomm, domain%topology, halo_width, 'full')
    call move_alloc(halo, halo_out)

end
end module ecs_halo_vec_c_factory_mod
