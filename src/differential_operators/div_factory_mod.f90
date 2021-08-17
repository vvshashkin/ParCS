module div_factory_mod

use partition_mod,     only : partition_t
use mesh_mod,          only : mesh_t
use abstract_div_mod,  only : div_operator_t
use div_c2_mod,        only : div_c2_t

implicit none

contains

function create_div_operator(mesh, partition, div_operator_name) result(div)
    type(partition_t), intent(in)    :: partition
    type(mesh_t),      intent(in)    :: mesh(partition%ts:partition%te)
    character(*),      intent(in)    :: div_operator_name

    class(div_operator_t), allocatable :: div

    if(div_operator_name == 'divergence_c2') then
        div = div_c2_t()
    elseif(div_operator_name == 'divergence_a2') then
        div = create_div_a2_operator(mesh,partition)
    else
        call avost("unknown divergence operator: "//div_operator_name)
    end if
end

function create_div_a2_operator(mesh, partition) result(div)

    use div_a2_mod,                 only : div_a2_t
    use mpi,                        only : MPI_comm_rank, MPI_comm_size, mpi_comm_world
    use exchange_factory_mod,       only : create_Agrid_halo_exchange
    use ecs_halo_mod,               only : ecs_halo_t
    use ecs_halo_factory_mod,       only : init_ecs_halo
    use ecs_halo_vec_a_factory_mod, only : init_ecs_halo_vect
    !use grid_function_factory_mod, only : create_grid_function

    type(partition_t), intent(in)    :: partition
    type(mesh_t),      intent(in)    :: mesh(partition%ts:partition%te)

    type(div_a2_t)  :: div

    class(ecs_halo_t), allocatable   :: halo
    integer(kind=4) :: myid, np, ierr, t

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    div%exch_halo = create_Agrid_halo_exchange( &
                      partition, 4, 'full', myid, np )

    allocate(div%halo_vec(partition%ts:partition%te))
    do t = partition%ts, partition%te
        halo = init_ecs_halo(mesh(t)%is, mesh(t)%ie, &
                             mesh(t)%js, mesh(t)%je, &
                             mesh(t)%nx, 2         , &
                             mesh(t)%hx)
        div%halo_vec(t) = init_ecs_halo_vect(mesh(t)%panel_ind,&
                                             mesh(t)%is, mesh(t)%ie,  &
                                             mesh(t)%js, mesh(t)%je,  &
                                             mesh(t)%nx, 1,    &
                                             mesh(t)%hx, halo)
     end do

end function create_div_a2_operator

end module div_factory_mod
