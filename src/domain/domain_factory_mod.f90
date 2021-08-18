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

    type(domain_t),   intent(out), target :: domain
    character(len=*), intent(in)          :: topology_type
    character(len=*), intent(in)          :: staggering_type
    integer(kind=4),  intent(in)          :: nh, nz

    !class(metric_t), allocatable  :: metric
    integer(kind=4) :: halo_width

!have to be passed as an argument in future
    halo_width = 8

    domain%topology = create_topology(topology_type)

    call create_metric(domain%metric,domain%topology,"ecs")

    call create_parcomm(parcomm_global%comm_w, domain%parcomm)

    call domain%partition%init(nh, nz, max(1,domain%parcomm%np/6), domain%parcomm%myid, domain%parcomm%Np, &
                                staggering_type, strategy = 'default')

    call create_mesh(domain%mesh_c, domain%partition, domain%metric, halo_width, 'c')
    call create_mesh(domain%mesh_x, domain%partition, domain%metric, halo_width, 'x')
    call create_mesh(domain%mesh_y, domain%partition, domain%metric, halo_width, 'y')
    call create_mesh(domain%mesh_xy, domain%partition, domain%metric, halo_width, 'xy')

    select case(staggering_type)
    case ('A')
        domain%mesh_p => domain%mesh_c
        domain%mesh_u => domain%mesh_c
        domain%mesh_v => domain%mesh_c
    case ('C')
        domain%mesh_p => domain%mesh_c
        domain%mesh_u => domain%mesh_x
        domain%mesh_v => domain%mesh_y
    case default
        call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//staggering_type)
    end select
end subroutine create_domain

end module domain_factory_mod
