module test_rk4_mod

implicit none

contains

subroutine test_rk4()
    use stvec_iomega_mod,    only: stvec_iomega_t, init_stvec_iomega
    use operator_iomega_mod, only: operator_iomega_t, init_operator_iomega
    use explicit_Eul1_mod,   only: explicit_Eul1_t
    use rk4_mod,             only: rk4_t, init_rk4
    use const_mod,           only : pi, Day24h_sec
    type(stvec_iomega_t) v1, v2
    type(operator_iomega_t) oper
    type(explicit_Eul1_t) ts_exEul
    type(rk4_t) ts_rk4
    integer, parameter :: N = 10
    real(kind=8), parameter :: dt = 6._8*3600._8
    complex(kind=8) omega(N), ftrue1(N), ftrue_rk4(N)
    real(kind=8) e1, e2
    real(kind=8), parameter :: tolerance = 1e-12_8

    integer i

    call init_stvec_iomega(v1,N)
    v1%f(1:N) = 1._8
    call v2%copy(v1)

    do i = 1,N
        omega(i) = cmplx(0._8, 2.0_8*pi*i/Day24h_sec)
        !Analytical responces of Eul1 & RK4 scheme 
        ftrue1(i) = 1._8+dt*omega(i)
        ftrue_rk4(i) = 1._8+dt*omega(i)+0.5_8*(dt*omega(i))**2      + &
                                              (dt*omega(i))**3/6._8 + &
                                              (dt*omega(i))**4/24._8
    end do

    call init_operator_iomega(oper, N, omega)

    call ts_exEul%step(oper, v1, dt)

    call init_rk4(ts_rk4, v2)
    call ts_rk4%step(oper, v2, dt)

    e1 = sum( abs( v1%f(1:N)-ftrue1(1:N)))
    e2 = sum( abs( v2%f(1:N)-ftrue_rk4(1:N)))
    if(e1 < tolerance) then
        print *, "explicit Euler 1st order scheme test passed"
    else
        print *, "explicit Euler 1st order scheme test failed"
    end if
    if(e2 < tolerance) then
        print *, "RK4 scheme test passed"
    else
        print *, "RK4 scheme test failed"
    end if

end subroutine test_rk4

end module test_rk4_mod
