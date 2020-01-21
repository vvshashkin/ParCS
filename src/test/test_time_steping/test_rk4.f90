module test_rk4_mod

implicit none

contains

subroutine test_rk4()
    use stvec_iomega_mod,    only: stvec_iomega_t, init_stvec_iomega
    use operator_iomega_mod, only: operator_iomega_t, init_operator_iomega
    use const_mod,           only : pi
    type(stvec_iomega_t) v1, v2
    type(operator_iomega_t) oper
    integer, parameter :: N = 10
    complex(kind=8) omega(N)

    integer i

    call init_stvec_iomega(v1,N)
    v1%f(1:N) = 1._8
    call v2%copy(v1)

    do i = 1,N
        omega(i) = (0._8, 2.0_8*pi*i/86400._8)
    end do

    call init_operator_iomega(oper, N, omega)

    call oper%act(v2,v1)

    print *, v1%f
    print *, "-------------"
    print *, v2%f

end subroutine test_rk4

end module test_rk4_mod
