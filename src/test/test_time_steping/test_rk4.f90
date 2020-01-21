module test_rk4_mod

implicit none

contains

subroutine test_rk4()
    use stvec_iomega_mod,    only: stvec_iomega_t, init_stvec_iomega
    use operator_iomega_mod, only: operator_iomega_t, init_operator_iomega
    use explicit_Eul1_mod,   only: explicit_Eul1_t
    use const_mod,           only : pi, Day24h_sec
    type(stvec_iomega_t) v1, v2
    type(operator_iomega_t) oper
    type(explicit_Eul1_t) ts_exEul
    integer, parameter :: N = 10
    real(kind=8), parameter :: dt = 6._8*3600._8
    complex(kind=8) omega(N)

    integer i

    call init_stvec_iomega(v1,N)
    v1%f(1:N) = 1._8
    call v2%copy(v1)

    do i = 1,N
        omega(i) = cmplx(0._8, 2.0_8*pi*i/Day24h_sec)
    end do

    call init_operator_iomega(oper, N, omega)

    !call oper%act(v2,v1)
    call ts_exEul%step(oper, v2, dt)

    print *, v1%f
    print *, "-------------"
    print *, v2%f

end subroutine test_rk4

end module test_rk4_mod
