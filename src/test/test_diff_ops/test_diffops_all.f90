program test_diffops

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, &
                                parcomm_global
use test_diffops_mod, only: err_container_t, test_div, test_grad, test_conv, test_curl, &
                            test_coriolis, test_KE, test_coriolis_vec_inv

implicit none

real(kind=8) :: err
type(err_container_t)  :: errs

call init_global_parallel_enviroment()

! errs = test_div(N=32,div_oper_name="divergence_a2_ecs",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_a2_cons",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_a2_fv",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_a2_fv"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_c2",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c2"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c_sbp21"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_c_sbp42"
!     print "(A,5E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ah2",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah2"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ah42_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah42_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_div(N=32,div_oper_name="divergence_ah43_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "divergence_ah43_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_ecs",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_a2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_a2_cons",staggering="A")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_a2_cons"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah21_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah21_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah42_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah42_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_ah43_sbp_ecs",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_ah43_sbp_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_c2_ecs",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c2_ecs"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! ! errs = test_grad(N=32,grad_oper_name="gradient_c2_cons",staggering="C")
! ! if(parcomm_global%myid == 0) then
! !     print *, "gradient_c2_cons"
! !     print "(A,4E15.7)", "Err: ", errs%values
! ! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_c_sbp21",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c_sbp21"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_grad(N=32,grad_oper_name="gradient_c_sbp42",staggering="C")
! if(parcomm_global%myid == 0) then
!     print *, "gradient_c_sbp42"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_curl(N=32,curl_oper_name="curl_divergence_ah42_sbp",staggering="Ah")
! if(parcomm_global%myid == 0) then
!     print *, "curl_divergence_ah42_sbp"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

errs = test_coriolis(N=32, coriolis_op_name="coriolis_colocated", staggering="A")
if (parcomm_global%myid==0) then
    print *, "coriolis_colocated"
    print "(A,4E15.7)", "Err: ", errs%values
end if

errs = test_coriolis_vec_inv(N=32, coriolis_op_name="coriolis_Cgrid", staggering="C")
if (parcomm_global%myid==0) then
    print *, "coriolis_Cgrid_vec_inv"
    print "(A,4E15.7)", "Err: ", errs%values
end if
!
! errs = test_KE(N=32, KE_oper_name="KE_colocated", staggering="A")
! if (parcomm_global%myid==0) then
!     print *, "KE_colocated"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if
!
! errs = test_KE(N=32, KE_oper_name="KE_Cgrid", staggering="C")
! if (parcomm_global%myid==0) then
!     print *, "KE_Cgrid"
!     print "(A,4E15.7)", "Err: ", errs%values
! end if

call deinit_global_parallel_enviroment()

end program test_diffops
