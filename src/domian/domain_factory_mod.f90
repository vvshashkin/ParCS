module domain_factory_mod

use domain_mod, only : domain_t
use mpi

implicit none

contains

subroutine create_ecs_global_domain(domain, nh, nz)

    use mesh_factory_mod, only : create_equiangular_mesh

    type(domain_t),  intent(out) :: domain
    integer(kind=4), intent(in)  :: nh, nz

    integer(kind=4) :: halo_width, myid, np, ts, te, ierr

    halo_width = 8

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    call domain%partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')

    ts = domain%partition%ts
    te = domain%partition%te

    call create_equiangular_mesh(domain%mesh_p, domain%partition, halo_width, 'C', 'p')
    call create_equiangular_mesh(domain%mesh_u, domain%partition, halo_width, 'C', 'u')
    call create_equiangular_mesh(domain%mesh_v, domain%partition, halo_width, 'C', 'v')

end subroutine create_ecs_global_domain

end module domain_factory_mod
