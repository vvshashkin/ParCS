module diag_NHlin_mod

use global_diag_mod, only : calc_global_diag!, integral_t

implicit none

    !type(integral_t) hmax, hmin, mass

contains

subroutine init_NHlin_diag_mod()

end subroutine init_NHlin_diag_mod

subroutine print_NHlin_diag(istep, stvec, model_params, myid, master_id)

    use mpi
    use stvec_NHlin_mod,        only : stvec_NHlin_t
    use parameters_NHlin_mod,   only : parameters_NHlin_t
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t

    integer(kind=4),          intent(in) :: istep
    type(stvec_NHlin_t),      intent(in) :: stvec
    type(parameters_NHlin_t), intent(in) :: model_params
    integer(kind=4),          intent(in) :: myid, master_id

    real(kind=8), allocatable :: tmax(:), tmin(:), mass(:)

    tmax = calc_global_diag(myid, master_id, .false., stvec, model_params, tmax_local, MPI_MAX)
    tmin = calc_global_diag(myid, master_id, .false., stvec, model_params, tmin_local, MPI_MIN)
    !mass = calc_global_diag(myid, master_id, .false., stvec, model_params, mass_local, MPI_SUM)

    if(myid==master_id) &
       print "(A,I7,A,2E15.7,A,E20.12)", "T=", istep, " THETA min/max=", tmin(1), tmax(1)!, &
                                               !" mass=", mass(1)

end subroutine print_NHlin_diag

function tmax_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
    use stvec_NHlin_mod,      only : stvec_NHlin_t
    use parameters_NHlin_mod, only : parameters_NHlin_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je, ks, ke

    select type(model_params)
    class is (parameters_NHlin_t)
    select type(stvec)
    class is (stvec_NHlin_t)

        allocate(f(1))
        ts = model_params%ts; te = model_params%te
        is  = model_params%tiles(ts)%is;    ie  = model_params%tiles(ts)%ie
        js  = model_params%tiles(ts)%js;    je  = model_params%tiles(ts)%je
        ks  = model_params%tiles(ts)%ks;    ke  = model_params%tiles(ts)%ke
        f(1) = maxval(stvec%theta(ts)%p(is:ie,js:je,ks-1:ke))

        do ind = ts+1, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            ks  = model_params%tiles(ind)%ks;    ke  = model_params%tiles(ind)%ke
            f(1) = max(maxval(stvec%theta(ind)%p(is:ie,js:je,ks-1:ke)), f(1))
        end do

    class default
        call avost("stvec type mismatch in hmax NHlin integral. stop!")
    end select
    class default
        call avost("model_params type mismatch in hmax NHlin integral. stop!")
    end select
end function tmax_local

function tmin_local(stvec, model_params) result(f)

    use mpi
    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
    use stvec_NHlin_mod,      only : stvec_NHlin_t
    use parameters_NHlin_mod, only : parameters_NHlin_t

    real(kind=8), allocatable                      :: f(:)
    class(state_abstract_t),            intent(in) :: stvec
    class(model_parameters_abstract_t), intent(in) :: model_params

    integer(kind=4) ts, te, ind, is, ie, js, je, ks, ke

    select type(model_params)
    class is (parameters_NHlin_t)
    select type(stvec)
    class is (stvec_NHlin_t)

        allocate(f(1))
        ts = model_params%ts; te = model_params%te
        is  = model_params%tiles(ts)%is;    ie  = model_params%tiles(ts)%ie
        js  = model_params%tiles(ts)%js;    je  = model_params%tiles(ts)%je
        ks  = model_params%tiles(ts)%ks;    ke  = model_params%tiles(ts)%ke
        f(1) = minval(stvec%theta(ts)%p(is:ie,js:je,ks-1:ke))

        do ind = ts+1, te
            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
            ks  = model_params%tiles(ind)%ks;    ke  = model_params%tiles(ind)%ke
            f(1) = min(minval(stvec%theta(ind)%p(is:ie,js:je,ks-1:ke)), f(1))
        end do

    class default
        call avost("stvec type mismatch in hmax NHlin integral. stop!")
    end select
    class default
        call avost("model_params type mismatch in hmax NHlin integral. stop!")
    end select
end function tmin_local

!function mass_local(stvec, model_params) result(f)
!
!    use mpi
!    use container_abstract_mod, only : state_abstract_t, model_parameters_abstract_t
!    use stvec_NHlin_mod,      only : stvec_NHlin_t
!    use parameters_NHlin_mod, only : parameters_NHlin_t
!
!    real(kind=8), allocatable                      :: f(:)
!    class(state_abstract_t),            intent(in) :: stvec
!    class(model_parameters_abstract_t), intent(in) :: model_params
!
!    integer(kind=4) ts, te, ind, is, ie, js, je, ks, ke, i, j, k
!
!    select type(model_params)
!    class is (parameters_NHlin_t)
!    select type(stvec)
!    class is (stvec_NHlin_t)
!
!        allocate(f(1))
!        ts = model_params%ts; te = model_params%te
!        f(1) = 0._8
!
!        do ind = ts, te
!            is  = model_params%tiles(ind)%is;    ie  = model_params%tiles(ind)%ie
!            js  = model_params%tiles(ind)%js;    je  = model_params%tiles(ind)%je
!            ks  = model_params%tiles(ind)%ks;    ke  = model_params%tiles(ind)%ke
!            do k = ks, ke
!                do j = js, je
!                    do i = is, ie
!                        f(1) = f(1) + stvec%h(ind)%p(i,j,k)*model_params%mesh(ind)%G(i,j)
!                    end do
!                end do
!            end do
!        end do
!
!    class default
!        call avost("stvec type mismatch in hmax NHlin integral. stop!")
!    end select
!    class default
!        call avost("model_params type mismatch in hmax NHlin integral. stop!")
!    end select
!end function mass_local

end module diag_NHlin_mod
