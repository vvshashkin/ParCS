module swlin_output_mod

use grid_function_mod, only : grid_function_t
use outputer_abstract_mod,       only : outputer_t

implicit none

type(grid_function_t), allocatable :: gf_buffer(:)
class(outputer_t),     allocatable :: outputer

contains

subroutine init_swlin_output(myid, master_id, np, partition, ts, te, mesh)
    use stvec_swlin_mod,      only : stvec_swlin_t
    use mesh_mod,             only : mesh_t
    use partition_mod,        only : partition_t
    use exchange_mod,         only : exchange_t
    use exchange_factory_mod, only : create_gather_exchange
    use outputer_factory_mod, only : create_master_paneled_outputer


    integer(kind=4),     intent(in) :: myid, master_id, np
    type(partition_t),   intent(in) :: partition
    integer(kind=4),     intent(in) :: ts, te
    type(mesh_t),        intent(in) :: mesh(ts:te)

    integer(kind=4) i
    type(exchange_t) exch_gather

    allocate(gf_buffer(ts:te))

    do i=ts, te
        call gf_buffer(i)%init(mesh(i)%panel_ind, mesh(i)%is, mesh(i)%ie,      &
                               mesh(i)%js, mesh(i)%je, mesh(i)%ks,mesh(i)%ke,  &
                               mesh(i)%halo_width,mesh(i)%halo_width, 0)
    end do

    call create_gather_exchange(exch_gather, partition, master_id, myid, np)
    outputer = create_master_paneled_outputer(master_id = master_id, gather_exch = exch_gather)
end subroutine init_swlin_output

subroutine write_swlin(myid, master_id, ts, te, stvec, ms, me, mesh, rec_num)
    use stvec_swlin_mod,      only : stvec_swlin_t
    use mesh_mod,             only : mesh_t

    integer(kind=4),     intent(in) :: myid, master_id
    integer(kind=4),     intent(in) :: ts, te
    type(stvec_swlin_t), intent(in) :: stvec
    integer(kind=4),     intent(in) :: ms, me
    type(mesh_t),        intent(in) :: mesh(ms:me)
    integer(kind=4),     intent(in) :: rec_num

    integer(kind=4) ind

    do ind = ts, te
        gf_buffer(ind) = stvec%h(ind)
    end do

    call outputer%write(gf_buffer, mesh, ms, me, "h.dat", rec_num)

    do ind = ts, te
        gf_buffer(ind) = stvec%u(ind)
    end do

    call outputer%write(gf_buffer, mesh, ms, me, "u.dat", rec_num)

    do ind = ts, te
        gf_buffer(ind) = stvec%v(ind)
    end do

    call outputer%write(gf_buffer, mesh, ms, me, "v.dat", rec_num)

end subroutine write_swlin

end module swlin_output_mod
