program latlon_output_test

    use test_latlon_output_mod,  only : test_latlon_output
    use parcomm_mod,             only : init_global_parallel_enviroment, &
                                        deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_latlon_output(staggering="A", scalar_grid="A")
    !call test_latlon_output(staggering="Ah",scalar_grid="Ah")
    !call test_latlon_output(staggering="C", scalar_grid="A")
    !call test_latlon_output(staggering="C", scalar_grid="Ah")

    call deinit_global_parallel_enviroment()

end program latlon_output_test
