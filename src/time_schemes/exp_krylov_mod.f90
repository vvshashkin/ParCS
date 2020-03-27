module exp_krylov_mod

use container_abstract_mod,  only : state_abstract_t, model_parameters_abstract_t
use stvec_abstract_mod,      only : stvec_abstract_t
use timescheme_abstract_mod, only : timescheme_abstract_t
use operator_abstract_mod,   only : operator_abstract_t

implicit none

type, extends(timescheme_abstract_t), public :: exp_krylov_t

    class(operator_abstract_t), allocatable :: operator
    class(stvec_abstract_t),    allocatable :: y(:)
    integer(kind=4)                         :: Mmax = 30 !Krylov subspace size
    integer(kind=4)                         :: iom  = 2 !Number of Krylov vectors to orthogonalize against

contains

    procedure, public :: step => step_exp_opt_collective

end type exp_krylov_t

contains

subroutine init_exp_krylov(ts_exp, operator, v, Mmax, iom)
    class(timescheme_abstract_t), allocatable, &
                                intent(out)   :: ts_exp
    class(operator_abstract_t), intent(in)    :: operator
    class(stvec_abstract_t),    intent(in)    :: v !example of model state vector
    integer(kind=4),  optional, intent(in)    :: Mmax, iom

    allocate(exp_krylov_t :: ts_exp)
    select type(ts_exp)
    class is (exp_krylov_t)
        allocate(ts_exp%y(1:ts_exp%Mmax+1), source=v)
        ts_exp%operator = operator
        if(present(Mmax)) ts_exp%Mmax = Mmax
        if(present(iom)) ts_exp%iom = iom
    end select

end subroutine init_exp_krylov

subroutine step_exp(this, v0, model_params, dt)

    use mpi

    class(exp_krylov_t),                intent(inout) :: this
    class(state_abstract_t),            intent(inout) :: v0
    class(model_parameters_abstract_t), intent(in)    :: model_params
    real(kind=8),                       intent(in)    :: dt

    integer(kind=4)  :: i, j
    real(8) h, beta
    real(8) Heis(this%Mmax,this%Mmax), ZZ1(this%Mmax, this%Mmax)

    select type(v0)
        class is (stvec_abstract_t)

    Heis = 0._8

    call this%y(1)%copy(v0)
    beta = this%y(1)%norm()
    call this%y(1)%smult(1._8/beta)

    do i=2,this%Mmax+1
        call this%operator%act(this%y(i),this%y(i-1),model_params)
        do j = max(1,i-this%iom),i-1
            h = this%y(i)%dot(this%y(j))
            call this%y(i)%add(this%y(j),1._8,-h)
            Heis(j,i-1) = h
        end do
        h = this%y(i)%norm()
        call this%y(i)%smult(1._8/h)
        if (i<=this%Mmax) Heis(i,i-1) = h
    end do

    call eig_exp(Heis,dt,this%Mmax,ZZ1)

    call v0%add(this%y(1),1._8,(ZZ1(1,1)-1._8)*beta)
    do i=2,this%Mmax
        call v0%add(this%y(i),1._8,ZZ1(i,1)*beta)
    end do

    class default
        call avost("EXP_Krylov scheme works only with class(stvec_abstract_t)")
    end select


end subroutine step_exp


subroutine step_exp_opt_collective(this, v0, model_params, dt)

    use mpi

    class(exp_krylov_t),                intent(inout) :: this
    class(state_abstract_t),            intent(inout) :: v0
    class(model_parameters_abstract_t), intent(in)    :: model_params
    real(kind=8),                       intent(in)    :: dt

    integer(kind=4)  :: i, j
    real(8) h(-this%iom:0), h_gl(-this%iom:0), beta, beta_loc
    real(8) Heis(this%Mmax,this%Mmax), ZZ1(this%Mmax, this%Mmax)
    integer(kind=4) ierr, myid

    select type(v0)
        class is (stvec_abstract_t)

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    Heis = 0._8
    h = 0._8
    h_gl = 0._8

    call this%y(1)%copy(v0)
    beta_loc = this%y(1)%norm()**2
    call mpi_allreduce(beta_loc,beta, 1, mpi_double, mpi_sum, mpi_comm_world, ierr)
    beta = sqrt(beta)
    call this%y(1)%smult(1._8/beta)

    do i=2,this%Mmax+1
        call this%operator%act(this%y(i),this%y(i-1),model_params)
        do j = max(1,i-this%iom),i-1
            h(j-i) = this%y(i)%dot(this%y(j))
        end do
        h(0) = this%y(i)%norm()**2

        call mpi_allreduce(h,h_gl, this%iom+1, mpi_double, mpi_sum, mpi_comm_world, ierr)
        h_gl(0) = sqrt(h_gl(0)-sum(h_gl(-this%iom:-1)**2))

        do j = max(1,i-this%iom),i-1
            Heis(j,i-1) = h_gl(j-i)
            call this%y(i)%add(this%y(j),1._8,-h_gl(j-i))
        end do
        call this%y(i)%smult(1._8/h_gl(0))
        if (i<=this%Mmax) Heis(i,i-1) = h_gl(0)
    end do

    call eig_exp(Heis,dt,this%Mmax,ZZ1)

    call v0%add(this%y(1),1._8,(ZZ1(1,1)-1._8)*beta)
    do i=2,this%Mmax
        call v0%add(this%y(i),1._8,ZZ1(i,1)*beta)
    end do

    class default
        call avost("EXP_Krylov scheme works only with class(stvec_abstract_t)")
    end select


end subroutine step_exp_opt_collective

subroutine eig_exp(H,dt,Mmax,expH)
    integer(kind=4) Mmax
    real(8) H(MMAX,MMAX)
    real(8) dt
    real(8) expH(MMAX,MMAX)
    real(8) P(MMAX,MMAX)
    real(8) WR(MMAX), WI(MMAX)
    real(8) FV(MMAX)
    complex(8) PP(MMAX,MMAX), P1(MMAX,MMAX), PR(Mmax,Mmax)
    complex(8) eval(Mmax)
    integer IA(MMAX)
    integer ierr, i, j
    complex(8), parameter :: im = (0._8,1._8)
    INTEGER IPVT(Mmax)
    real(8) DET(2)
    complex(8) ZWK(2*Mmax), ZDET(2)
    real(8) Rcond
    real(8), parameter :: zeps = 1e-14_8
    real(8) z

    expH = H
    call rg(Mmax,Mmax,expH,WR,WI,1,P,IA,FV,ierr)

    j = 1

    !print *, "eig", maxval(wr), minval(wr), maxval(wi), minval(wi)
    do while(j<=Mmax)
        if(abs(wi(j))<1d-15) then
            PP(1:Mmax,j) = P(1:Mmax,j)
            eval(j) = wr(j)+im*wi(j)
            j = j+1
        else
            PP(1:Mmax,j)   = P(1:Mmax,j)+im*P(1:Mmax,j+1)
            PP(1:Mmax,j+1) = P(1:Mmax,j)-im*P(1:Mmax,j+1)
            eval(j) = wr(j)+im*wi(j)
            eval(j+1) = wr(j)-im*wi(j)
            j = j+2
        end if
    end do
    do j=1, Mmax
        if(abs(eval(j))> zeps) then
            !eval(j) = (exp(eval(j)*dt)-1._8)/(eval(j)*dt)
            eval(j) = exp(eval(j)*dt)
        else
            eval(j) = 1._8
        end if
    end do

    P1 = PP
    CALL ZGECO(P1,Mmax,Mmax,IPVT,RCOND,ZWK)
    CALL ZGEDI(P1,Mmax,Mmax,IPVT,ZDET,ZWK,11)

    do j=1,Mmax
      do i=1,Mmax
        expH(i,j) = real(sum(PP(i,:)*eval(:)*P1(:,j)))
      end do
    end do

end subroutine eig_exp

end module exp_krylov_mod
