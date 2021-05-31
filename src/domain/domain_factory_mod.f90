module domain_factory_mod

use domain_mod, only : domain_t
use topology_factory_mod,      only: init_topology
use cubed_sphere_topology_mod, only: cubed_sphere_topology_t
use mpi

implicit none

contains

subroutine create_domain(domain, topology_type, staggering_type, nh, nz)

    use mesh_factory_mod,    only : create_equiangular_mesh
    use parcomm_factory_mod, only : create_parcomm

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz

    integer(kind=4) :: halo_width

!have to be passed as an argument in future
    halo_width = 8

    domain%topology = init_topology(topology_type)

    call create_parcomm(domain%parcomm)

    call domain%partition%init(nh, nz, max(1,domain%parcomm%np/6), domain%parcomm%myid, domain%parcomm%Np, &
                                staggering_type, strategy = 'default')

    call create_equiangular_mesh(domain%mesh_p, domain%partition, halo_width, staggering_type, 'p')
    call create_equiangular_mesh(domain%mesh_u, domain%partition, halo_width, staggering_type, 'u')
    call create_equiangular_mesh(domain%mesh_v, domain%partition, halo_width, staggering_type, 'v')

end subroutine create_domain

end module domain_factory_mod
