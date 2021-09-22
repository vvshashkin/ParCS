program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment
use test_diffops_mod, only: err_container_t, test_div, test_grad, test_conv, test_curl

implicit none

real(kind=8) :: err
type(err_container_t)  :: errs
integer(kind=4), parameter :: Ns(3) = [32,64,128]

call init_global_parallel_enviroment()

! errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
! print *, "divergence_a2_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_a2_cons",staggering="A")
! print *, "divergence_a2_cons"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_a2_fv",staggering="A")
! print *, "divergence_a2_fv"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_c2",staggering="C")
! print *, "divergence_c2"
! print "(A,5E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_c_sbp21",staggering="C")
! print *, "divergence_c_sbp21"
! print "(A,5E15.7)", "Err: ", errs%values
!
 errs = test_div(N=32,div_oper_name="divergence_c_sbp42",staggering="C")
 print *, "divergence_c_sbp42"
 print "(A,5E15.7)", "Err: ", errs%values
!
!
! errs = test_div(N=32,div_oper_name="divergence_ah2",staggering="Ah")
! print *, "divergence_ah2"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_ah42_sbp",staggering="Ah")
! print *, "divergence_ah42_sbp"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_div(N=32,div_oper_name="divergence_ah43_sbp",staggering="Ah")
! print *, "divergence_ah43_sbp"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_ecs",staggering="A")
! print *, "gradient_a2_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_cons",staggering="A")
! print *, "gradient_a2_cons"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah2_ecs",staggering="Ah")
! print *, "gradient_ah2_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=64,grad_oper_name="gradient_ah42_sbp_ecs",staggering="Ah")
! print *, "gradient_ah42_sbp_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=64,grad_oper_name="gradient_ah43_sbp_ecs",staggering="Ah")
! print *, "gradient_ah43_sbp_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_c2_ecs",staggering="C")
! print *, "gradient_c2_ecs"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_c2_cons",staggering="C")
! print *, "gradient_c2_cons"
! print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_grad(N=32,grad_oper_name="gradient_c_sbp21",staggering="C")
! print *, "gradient_c_sbp21"
! print "(A,4E15.7)", "Err: ", errs%values
!
 errs = test_grad(N=32,grad_oper_name="gradient_c_sbp42",staggering="C")
 print *, "gradient_c_sbp42"
 print "(A,4E15.7)", "Err: ", errs%values
!
! errs = test_curl(N=32,div_oper_name="divergence_ah42_sbp",staggering="Ah")
! print *, "curl_divergence_ah42_sbp"
! print "(A,4E15.7)", "Err: ", errs%values

! call test_conv(operator_name="gradient_c_sbp21",staggering="C",Ns=Ns)
! call test_conv(operator_name="gradient_c_sbp42",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp21",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_c_sbp42",staggering="C",Ns=Ns)
! call test_conv(operator_name="divergence_ah42_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="divergence_ah43_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah42_sbp",staggering="Ah",Ns=Ns)
! call test_conv(operator_name="curl_divergence_ah43_sbp",staggering="Ah",Ns=Ns)

call deinit_global_parallel_enviroment()

end program test_diffops
