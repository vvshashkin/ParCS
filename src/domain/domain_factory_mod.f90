module domain_factory_mod

use domain_mod,                only : domain_t
use topology_factory_mod,      only : create_topology
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use metric_mod,                only : metric_t
use metric_factory_mod,        only : create_metric
use mpi

implicit none

contains

subroutine create_domain(domain, topology_type, staggering_type, nh, nz)

    use mesh_factory_mod,    only : create_mesh
    use parcomm_factory_mod, only : create_parcomm
    use parcomm_mod,         only : parcomm_global

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz

    class(metric_t), allocatable  :: metric
    integer(kind=4) :: halo_width

!have to be passed as an argument in future
    halo_width = 8

    domain%topology = create_topology(topology_type)

    call create_metric(metric,domain%topology,"ecs")

    call create_parcomm(parcomm_global%comm_w, domain%parcomm)

    call domain%partition%init(nh, nz, max(1,domain%parcomm%np/6), domain%parcomm%myid, domain%parcomm%Np, &
                                staggering_type, strategy = 'default')

    call create_mesh(domain%mesh_p, domain%partition, metric, halo_width, staggering_type, 'p')
    call create_mesh(domain%mesh_u, domain%partition, metric, halo_width, staggering_type, 'u')
    call create_mesh(domain%mesh_v, domain%partition, metric, halo_width, staggering_type, 'v')

end subroutine create_domain

end module domain_factory_mod
