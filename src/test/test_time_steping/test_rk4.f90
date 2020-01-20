module test_rk4_mod

implicit none

contains

subroutine test_rk4()
    use stvec_iomega_mod, only: stvec_iomega_t, init_stvec_iomega
    type(stvec_iomega_t) v1, v2
    integer, parameter :: N = 10

    call init_stvec_iomega(v1,N)
    v1%f(1:N) = 1._8
    call v2.copy(v1)
    v1%f(1:N) = (1._8,2._8)
    call v1%add(v2,2._8,-2._8)

    print *, v1%f

end subroutine test_rk4

end module test_rk4_mod
