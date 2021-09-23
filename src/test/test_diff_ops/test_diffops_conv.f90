program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, &
                                parcomm_global
use test_diffops_mod, only: err_container_t, test_div, test_grad, test_conv, test_curl

implicit none

real(kind=8) :: err
type(err_container_t)  :: errs
integer(kind=4), parameter :: Ns(3) = [32,64,128]

call init_global_parallel_enviroment()

call test_conv(operator_name="gradient_c_sbp21",staggering="C",Ns=Ns)
call test_conv(operator_name="gradient_c_sbp42",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp21",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp42",staggering="C",Ns=Ns)
call test_conv(operator_name="divergence_ah42_sbp",staggering="Ah",Ns=Ns)
call test_conv(operator_name="divergence_ah43_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah42_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah43_sbp",staggering="Ah",Ns=Ns)

call deinit_global_parallel_enviroment()

end program test_diffops
