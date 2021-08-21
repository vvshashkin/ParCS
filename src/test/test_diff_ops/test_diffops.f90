program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment
use test_diffops_mod, only: err_container_t, test_div, test_grad_a2

implicit none

real(kind=8) :: err
type(err_container_t)  :: errs

call init_global_parallel_enviroment()

!errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
!print "(A,4E15.7)", "Err: ", errs%values

errs = test_div(N=32,div_oper_name="divergence_ah2",staggering="Ah")
print "(A,4E15.7)", "Err", errs%values

!err = test_grad_a2(N=64)
!print *, "Err: ", err

call deinit_global_parallel_enviroment()

end program test_diffops
