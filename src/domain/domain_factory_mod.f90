module domain_factory_mod

use domain_mod,                only : domain_t
use topology_factory_mod,      only : create_topology
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use metric_mod,                only : metric_t
use metric_factory_mod,        only : create_metric
use mpi

implicit none

contains

subroutine create_domain(domain, topology_type, staggering_type, nh, nz, &
                         parcomm)

    use mesh_factory_mod,    only : create_mesh
    use parcomm_factory_mod, only : create_parcomm
    use parcomm_mod,         only : parcomm_global, parcomm_t

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz
    type(parcomm_t),  optional, intent(in) :: parcomm


    !class(metric_t), allocatable  :: metric
    integer(kind=4) :: halo_width
    integer(kind=4) :: mpi_comm_local

!have to be passed as an argument in future
    halo_width = 8

    domain%topology = create_topology(topology_type)

    call create_metric(domain%metric,domain%topology,"ecs")

    domain%horizontal_staggering = staggering_type

    mpi_comm_local = parcomm_global%comm_w
    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        !call create_parcomm(parcomm_glocal%comm_w, domain%parcomm)
        domain%parcomm = parcomm_global
    end if

    call domain%partition%init(nh, nz, max(1,domain%parcomm%np/6), domain%parcomm%myid, domain%parcomm%Np, &
                               staggering_type, strategy = 'default')

    call create_mesh(domain%mesh_o, domain%partition, domain%metric, halo_width, 'c')
    call create_mesh(domain%mesh_x, domain%partition, domain%metric, halo_width, 'x')
    call create_mesh(domain%mesh_y, domain%partition, domain%metric, halo_width, 'y')
    call create_mesh(domain%mesh_xy, domain%partition, domain%metric, halo_width, 'xy')

    select case(staggering_type)
    case ('A')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_o
        domain%mesh_v = domain%mesh_o
    case ('Ah') !all degrees of freedom at corner points
        domain%mesh_p = domain%mesh_xy
        domain%mesh_u = domain%mesh_xy
        domain%mesh_v = domain%mesh_xy
    case ('C')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_x
        domain%mesh_v = domain%mesh_y
    case default
        call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//staggering_type)
    end select
end subroutine create_domain

end module domain_factory_mod
