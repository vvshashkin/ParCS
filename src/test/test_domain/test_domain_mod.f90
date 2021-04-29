module test_domain_mod

use domain_mod, only : domain_t


implicit none

contains

subroutine test_domain()

    use domain_factory_mod, only : create_ecs_global_domain

    type(domain_t) :: domain

    integer(kind=4) :: nh=100, nz=10, halo_width=50

    call create_ecs_global_domain(domain, nh, nz)


end subroutine test_domain

end module test_domain_mod
