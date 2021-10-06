program test_quadrature_main

    use parcomm_mod,         only : init_global_parallel_enviroment, &
                                    deinit_global_parallel_enviroment

    use test_quadrature_mod, only : test_quadrature

    implicit none

    call init_global_parallel_enviroment()

    call test_quadrature("default_quadrature", max_order_exact = 1, staggering="A", points_type="o")

    call deinit_global_parallel_enviroment()


end program test_quadrature_main
