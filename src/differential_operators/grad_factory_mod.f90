module grad_factory_mod

use abstract_grad_mod, only : grad_operator_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_grad_operator(domain, grad_operator_name) result(grad)
    type(domain_t),    intent(in)  :: domain
    character(len=*),  intent(in)  :: grad_operator_name

    class(grad_operator_t), allocatable :: grad

    if(grad_operator_name == 'gradient_c2') then
       ! grad = create_grad_contra_c2_operator(mesh,partition,grad_operator_name)
        call parcomm_global%abort("not implemented "//grad_operator_name)
    else if(grad_operator_name == 'gradient_a2_ecs' .or. &
            grad_operator_name == 'gradient_a2_cons') then
        grad = create_grad_contra_a2_operator(domain,grad_operator_name)
    else
        call parcomm_global%abort("unknown gradient operator: "//grad_operator_name)
    end if
end

! function create_grad_contra_c2_operator(mesh, partition,grad_operator_name) result(grad)

!     use grad_contra_c2_mod,        only : grad_contra_c2_t
!     use mpi,                       only : MPI_comm_rank, MPI_comm_size, mpi_comm_world
!     use exchange_factory_mod,      only : create_Agrid_halo_exchange
!     use ecs_halo_factory_mod,      only : init_ecs_halo
!     use ecs_halo_lin_factory_mod,  only : init_ecs_halo_lin
!     !use grid_function_factory_mod, only : create_grid_function

!     type(partition_t), intent(in)    :: partition
!     type(mesh_t),      intent(in)    :: mesh(partition%ts:partition%te)
!     character(*),      intent(in)    :: grad_operator_name

!     type(grad_contra_c2_t) :: grad

!     integer(kind=4) :: myid, np, ierr, t

!     call MPI_comm_rank(mpi_comm_world , myid, ierr)
!     call MPI_comm_size(mpi_comm_world , Np  , ierr)

!     grad%kernel_type = grad_operator_name

!     grad%exch_halo = create_Agrid_halo_exchange( &
!                       partition, 4, 'full', myid, np )

!     allocate(grad%halo(partition%ts:partition%te))
!     do t = partition%ts, partition%te
!          grad%halo(t) = init_ecs_halo_lin(mesh(t)%is, mesh(t)%ie, &
!                                       mesh(t)%js, mesh(t)%je, &
!                                       mesh(t)%nx, 2         , &
!                                       mesh(t)%hx)
!     end do

! end function create_grad_contra_c2_operator

function create_grad_contra_a2_operator(domain, grad_operator_name) result(grad)

    use grad_contra_a2_mod, only : grad_contra_a2_t
    use halo_factory_mod,   only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_contra_a2_t) :: grad

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    grad = grad_contra_a2_t()
    if(grad_operator_name=="gradient_a2_ecs") then
        call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")
    else if(grad_operator_name=="gradient_a2_cons") then
        call create_halo_procedure(grad%halo_procedure,domain,default_halo_width,"A_default")
    else
        call parcomm_global%abort("grad_factory_mod, create_grad_contra_a2_operator "//&
                                  "unknown gradient_a2 subtype: "// grad_operator_name)
    end if

end function create_grad_contra_a2_operator

end module grad_factory_mod
