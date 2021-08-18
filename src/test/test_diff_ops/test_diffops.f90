program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment
use test_diffops_mod, only: test_div_a2, test_grad_a2, test_laplace_spectre

implicit none

real(kind=8) :: err

call init_global_parallel_enviroment()

err = test_div_a2(N=32)
print *, "Err: ", err

err = test_grad_a2(N=64)
print *, "Err: ", err

call deinit_global_parallel_enviroment()

end program test_diffops