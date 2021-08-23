program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment
use test_diffops_mod, only: err_container_t, test_div, test_grad

implicit none

real(kind=8) :: err
type(err_container_t)  :: errs

call init_global_parallel_enviroment()

!errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
!print "(A,4E15.7)", "Err: ", errs%values

!errs = test_div(N=64,div_oper_name="divergence_ah2",staggering="Ah")
!print "(A,4E15.7)", "Err", errs%values

errs = test_grad(N=64,grad_oper_name="gradient_ah2_ecs",staggering="Ah")
print "(A,4E15.7)", "Err: ", errs%values

call deinit_global_parallel_enviroment()

end program test_diffops
