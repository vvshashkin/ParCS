module test_gl_diag_mod

use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
use partition_mod,          only : partition_t
use grid_function_mod,      only : grid_function_t
use tile_mod,               only : tile_t

implicit none

type, extends(model_parameters_abstract_t) :: diag_test_param_t
    integer(kind=4)                    :: nx, nz, npanels, halo_width
    integer(kind=4)                    :: ts, te
    type(partition_t)                  :: partition
    type(tile_t), allocatable          :: tiles(:)
end type diag_test_param_t

type, extends(state_abstract_t) :: diag_test_state_t
    type(grid_function_t), allocatable :: gf(:)
end type diag_test_state_t

contains

subroutine test_gl_diag()
    use mpi
    use parameters_swlin_mod,    only : parameters_swlin_t, init_swlin_parameters
    use stvec_swlin_mod,         only : stvec_swlin_t, init_stvec_swlin
    use global_diag_mod,         only : calc_global_diag!, integral_t


    type(diag_test_param_t)                   :: params
    type(diag_test_state_t)                   :: state
    integer(kind=4)                           :: myid, Np, ierr

    integer(kind=4), parameter                :: master_id = 0
    integer(kind=4), parameter                :: nx = 32, nz = 2
    integer(kind=4), parameter                :: halo_width = 0, npanels = 6

    integer(kind=4)                           :: is, ie, js, je
    integer(kind=4)                           :: ks, ke, panel_ind
    integer(kind=4)                           :: ind, i, j, k
    real(kind=8), allocatable                 :: ave_n_disp(:)

    logical lpass
    real(kind=8), parameter :: tolerance = 1e-15_8
    real(kind=8) true_answer(2*nz)

    call MPI_comm_rank(mpi_comm_world , myid, ierr)
    call MPI_comm_size(mpi_comm_world , Np  , ierr)

    !Init parameters
    params%nx = nx; params%nz = nz; params%npanels = npanels
    call params%partition%init(params%nx, params%nz, max(1,Np/6), myid, Np, strategy = 'default')
    params%ts = findloc(params%partition%proc_map, myid, dim=1)
    params%te = findloc(params%partition%proc_map, myid, back = .true., dim=1)
    allocate(params%tiles(params%ts:params%te))
    params%tiles(params%ts:params%te) = params%partition%tile(params%ts:params%te)

    !Init state
    allocate(state%gf(params%ts:params%te))
    do ind = params%ts, params%te

        panel_ind = params%tiles(ind)%panel_number
        is = params%tiles(ind)%is;     ie = params%tiles(ind)%ie
        js = params%tiles(ind)%js;     je = params%tiles(ind)%je
        ks = params%tiles(ind)%ks;     ke = params%tiles(ind)%ke

        call state%gf(ind)%init(panel_ind,is, ie, js, je, ks, ke, &
                                halo_width, halo_width, 0)
        do k=ks, ke
            do j=js, je
                do i = is, ie
                    state%gf(ind)%p(i,j,k) = (-1._8)**(k*i)
                end do
            end do
        end do

    end do

    do k=1,nz
        if(mod(k,2) == 0) then
            true_answer(2*k-1:2*k) = [1._8, 0._8]
        else
            true_answer(2*k-1:2*k) = [0._8, 1._8]
        end if
    end do

    ave_n_disp = calc_global_diag(myid, master_id, .true., state, params, &
                                  sum_local, MPI_SUM, get_ave_n_disp)
    lpass = maxval(abs(true_answer(1:2*nz)-ave_n_disp(1:2*nz))) < tolerance
    if(lpass) then
        print *, "global diag all reduce test passed on proc n=", myid
    else
        print *, "global diag all reduce test failed on proc n=", myid
    end if

    deallocate(ave_n_disp)
    ave_n_disp = calc_global_diag(myid, master_id, .false., state, params, &
                                  sum_local, MPI_SUM, get_ave_n_disp)
    if(myid == master_id) then
        lpass = maxval(abs(true_answer(1:2*nz)-ave_n_disp(1:2*nz))) < tolerance
        if(lpass) then
            print *, "global diag master reduce test passed"
        else
            print *, "global diag master reduce test failed"
        end if
    end if

end subroutine test_gl_diag

function sum_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je, i, j, nz, k
    real(kind=8) zfn

    select type(model_params)
    class is (diag_test_param_t)
    select type(stvec)
    class is (diag_test_state_t)

        zfn = 1._8 / (model_params%nx**2*model_params%npanels)
        nz = model_params%nz
        allocate(f(2*nz))
        ts = model_params%ts; te = model_params%te
        f(1:2*nz) = 0._8

        do ind = ts, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            do k=1, nz
                do j = js, je
                    do i = is, ie
                        f(2*k-1) = f(2*k-1) + stvec%gf(ind)%p(i,j,k)
                        f(2*k)   = f(2*k)   + stvec%gf(ind)%p(i,j,k)**2
                    end do
                end do
            end do
        end do

        f(1:2*nz) = f(1:2*nz)*zfn

    class default
        call avost("diag test failed: state type mismatch. stop!")
    end select
    class default
        call avost("diag test failed: param type mismatch. stop!")
    end select
end function sum_local

function get_ave_n_disp(f, model_params) result(q)
    real(kind=8), allocatable                      :: q(:)
    real(kind=8), intent(in)                       :: f(:)
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer nz2, k
    nz2 = size(f)
    allocate(q(1:nz2))
    q(1:nz2) = f(1:nz2)

    do k=2, nz2, 2
        q(k) = q(k)-q(k-1)**2 !no sqrt to avoid round-off errors when disp = sqrt(1-1)
    end do

end function get_ave_n_disp

end module test_gl_diag_mod
