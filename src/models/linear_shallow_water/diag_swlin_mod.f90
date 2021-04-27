module diag_swlin_mod

use global_diag_mod, only : calc_global_diag!, integral_t

implicit none

    !type(integral_t) hmax, hmin, mass

contains

subroutine init_swlin_diag_mod()

end subroutine init_swlin_diag_mod

subroutine print_swlin_diag(istep, stvec, model_params, myid, master_id)

    use mpi
    use stvec_swlin_mod,        only : stvec_swlin_t
    use parameters_swlin_mod,   only : parameters_swlin_t
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t

    integer(kind=4),          intent(in) :: istep
    type(stvec_swlin_t),      intent(in) :: stvec
    type(parameters_swlin_t), intent(in) :: model_params
    integer(kind=4),          intent(in) :: myid, master_id

    real(kind=8), allocatable :: hmax(:), hmin(:), mass(:)

    hmax = calc_global_diag(myid, master_id, .false., stvec, model_params, hmax_local, MPI_MAX)
    hmin = calc_global_diag(myid, master_id, .false., stvec, model_params, hmin_local, MPI_MIN)
    mass = calc_global_diag(myid, master_id, .false., stvec, model_params, mass_local, MPI_SUM)

    if(myid==master_id) &
       print "(A,I7,A,2E15.7,A,E20.12)", "T=", istep, " H min/max=", hmin(1), hmax(1), &
                                               " mass=", mass(1)

end subroutine print_swlin_diag

function hmax_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
    use stvec_swlin_mod,      only : stvec_swlin_t
    use parameters_swlin_mod, only : parameters_swlin_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je

    select type(model_params)
    class is (parameters_swlin_t)
    select type(stvec)
    class is (stvec_swlin_t)

        allocate(f(1))
        ts = model_params%ts; te = model_params%te
        is  = model_params%tiles(ts)%is;    ie  = model_params%tiles(ts)%ie
        js  = model_params%tiles(ts)%js;    je  = model_params%tiles(ts)%je
        f(1) = maxval(stvec%h%block(ts)%p(is:ie,js:je,1))

        do ind = ts+1, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            f(1) = max(maxval(stvec%h%block(ind)%p(is:ie,js:je,1)), f(1))
        end do

    class default
        call avost("stvec type mismatch in hmax swlin integral. stop!")
    end select
    class default
        call avost("model_params type mismatch in hmax swlin integral. stop!")
    end select
end function hmax_local

function hmin_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
    use stvec_swlin_mod,      only : stvec_swlin_t
    use parameters_swlin_mod, only : parameters_swlin_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je

    select type(model_params)
    class is (parameters_swlin_t)
    select type(stvec)
    class is (stvec_swlin_t)

        allocate(f(1))
        ts = model_params%ts; te = model_params%te
        is  = model_params%tiles(ts)%is;    ie  = model_params%tiles(ts)%ie
        js  = model_params%tiles(ts)%js;    je  = model_params%tiles(ts)%je
        f(1) = minval(stvec%h%block(ts)%p(is:ie,js:je,1))

        do ind = ts+1, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            f(1) = min(minval(stvec%h%block(ind)%p(is:ie,js:je,1)), f(1))
        end do

    class default
        call avost("stvec type mismatch in hmax swlin integral. stop!")
    end select
    class default
        call avost("model_params type mismatch in hmax swlin integral. stop!")
    end select
end function hmin_local

function mass_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
    use stvec_swlin_mod,      only : stvec_swlin_t
    use parameters_swlin_mod, only : parameters_swlin_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je, i, j

    select type(model_params)
    class is (parameters_swlin_t)
    select type(stvec)
    class is (stvec_swlin_t)

        allocate(f(1))
        ts = model_params%ts; te = model_params%te
        f(1) = 0._8

        do ind = ts, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            do j = js, je
                do i = is, ie
                    f(1) = f(1) + stvec%h%block(ind)%p(i,j,1)*model_params%mesh(ind)%G(i,j)
                end do
            end do
        end do

    class default
        call avost("stvec type mismatch in hmax swlin integral. stop!")
    end select
    class default
        call avost("model_params type mismatch in hmax swlin integral. stop!")
    end select
end function mass_local

end module diag_swlin_mod
