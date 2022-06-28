program test_diffops_3d

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_diffops_3d_mod, only : test_w2uv_interp, test_uv2w_interp, test_scalar_advection_3d, &
                                test_grad_3d, test_div_3d, test_co2contra_3d
use key_value_mod,       only : key_value_r8_t

implicit none

type(key_value_r8_t)  :: errs
character(len=:), allocatable :: config_str

call init_global_parallel_enviroment()

errs =  test_div_3d(Nh = 32, Nz = 20, &
                    hor_div_name = "divergence_c_sbp42", &
                    diff_eta_name = "eta_diff_w2p_sbp42", &
                    horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "div_3d_c_sbp42"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_div_3d(Nh = 32, Nz = 20, &
                    hor_div_name = "divergence_ah42_sbp", &
                    diff_eta_name = "eta_diff_w2p_sbp21", &
                    horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "div_3d_Ah_sbp42_21"
    print "(A,4E25.16)", "Err: ", errs%values
end if

! errs =  test_grad_3d(Nh = 64, Nz = 20, &
!                     hor_grad_name = "gradient_c_sbp42", &
!                     diff_eta_name = "eta_diff_p2w_sbp42", &
!                     horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "grad_3d_c_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_grad_3d(Nh = 32, Nz = 8, &
!                     hor_grad_name = "gradient_ah42_sbp_ecs", &
!                     diff_eta_name = "eta_diff_p2w_sbp21", &
!                     horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "grad_3d_Ah_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_colocated", &
!                     horizontal_staggering = "Ah", vertical_staggering = "None")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_colocated"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_Cgrid_h_sbp42_v_sbp42", &
!                     horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_C_CharneyPhilips"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
! errs =  test_co2contra_3d(Nh = 32, Nz = 8, &
!                     co2contra_3d_oper_name = "co2contra_3d_h_colocated_v_sbp42", &
!                     horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "co2contra_3d_h_colocated_v_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_colocated", &
!                          w2uv_hor_part_name     = "None",           &
!                          w2uv_vert_part_name    = "None",           &
!                          horizontal_staggering  = "Ah", vertical_staggering = "None")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_colocated"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_hor_colocated",         &
!                          w2uv_hor_part_name     = "None",                       &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp21",  &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_hor_colocated_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_hor_colocated",         &
!                          w2uv_hor_part_name     = "None",                       &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp42",  &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_hor_colocated_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_staggered",             &
!                          w2uv_hor_part_name     = "interp2d_p2uv_C_sbp42",     &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp21",  &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_staggered_C_sbp42_v_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_w2uv_interp(Nh = 32, Nz = 8, &
!                          w2uv_interpolator_name = "w2uv_staggered",             &
!                          w2uv_hor_part_name     = "interp2d_p2uv_C_sbp42",      &
!                          w2uv_vert_part_name    = "vertical_interp_w2p_sbp42",  &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "w2uv_staggered_C_sbp42_v_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
!                          uv2w_interpolator_name = "uv2w_colocated", &
!                          uv2w_hor_part_name     = "",               &
!                          uv2w_vert_part_name    = "",               &
!                          horizontal_staggering = "Ah", vertical_staggering = "None")
! if (parcomm_global%myid==0) then
!     print *, "uv2w_colocated"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
!                          uv2w_interpolator_name = "uv2w_hor_colocated",  &
!                          uv2w_hor_part_name     = "",                          &
!                          uv2w_vert_part_name    = "vertical_interp_p2w_sbp21", &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "uv2w_hor_colocated_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
!                          uv2w_interpolator_name = "uv2w_hor_colocated",  &
!                          uv2w_hor_part_name     = "",                          &
!                          uv2w_vert_part_name    = "vertical_interp_p2w_sbp42", &
!                          horizontal_staggering = "Ah", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "uv2w_hor_colocated_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
!                          uv2w_interpolator_name = "uv2w_staggered",                  &
!                          uv2w_hor_part_name     = "interp2d_uv2pvec_C_sbp42",           &
!                          uv2w_vert_part_name    = "vertical_interp_p2w_sbp21",       &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "uv2w_staggered_C_sbp42_v_sbp21"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_uv2w_interp(Nh = 32, Nz = 8, &
!                          uv2w_interpolator_name = "uv2w_staggered",                  &
!                          uv2w_hor_part_name     = "interp2d_uv2pvec_C_sbp42",           &
!                          uv2w_vert_part_name    = "vertical_interp_p2w_sbp42",       &
!                          horizontal_staggering = "C", vertical_staggering = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "uv2w_staggered_C_sbp42_v_sbp42"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
!                                  advection_oper_name      = "advection_w_staggered", &
!                                  config_str               = "&w_advection_conf "          // &
!                                  "hor_advection_oper_name = 'up4',"                       // &
!                                  "z_advection_oper_name   = 'adv_z_c2',"                  // &
!                                  "uv2w_operator_name      = 'uv2w_staggered',"            // &
!                                  "uv2w_hor_part_name      = 'interp2d_uv2pvec_C_sbp42'  " // &
!                                  "uv2w_vert_part_name     = 'vertical_interp_p2w_sbp42' " // &
!                                  "w_halo                  = 'ECS_Oz' /",                     &
!                                  points_type              = "w",                     &
!                                  horizontal_staggering    = "C",                     &
!                                  vertical_staggering      = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "advection_w_staggered, up4, adv_z_c2"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

! errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
!                                  advection_oper_name      = "advection_w_Ah", &
!                                  config_str               = "&w_advection_conf "          // &
!                                  "hor_advection_oper_name = 'sbp_d42',"                   // &
!                                  "z_advection_oper_name   = 'adv_z_c2',"                  // &
!                                  "uv2w_operator_name      = 'uv2w_hor_colocated',"       // &
!                                  "uv2w_vert_part_name     = 'vertical_interp_p2w_sbp42' /",                     &
!                                  points_type              = "w",                      &
!                                  horizontal_staggering    = "Ah",                     &
!                                  vertical_staggering      = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "advection_w_Ah, sbp_d42, adv_z_c2"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if

! errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
!                                  advection_oper_name      = "advection_p_staggered", &
!                                  config_str               = "&p_advection_conf "          // &
!                                  "hor_advection_oper_name = 'up4',"                       // &
!                                  "z_advection_oper_name   = 'adv_z_c2',"                  // &
!                                  "w2p_operator_name       = 'vertical_interp_w2p_sbp42'," // &
!                                  "uv2p_operator_name      = 'interp2d_uv2pvec_C_sbp42'  " // &
!                                  "p_halo                  = 'ECS_O' /",                      &
!                                  points_type              = "p",                     &
!                                  horizontal_staggering    = "C",                     &
!                                  vertical_staggering      = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "advection_p_staggered, up4, adv_z_c2"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if
!
! errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
!                                  advection_oper_name      = "advection_p_Ah", &
!                                  config_str               = "&p_advection_conf "          // &
!                                  "hor_advection_oper_name = 'sbp_d42',"                   // &
!                                  "z_advection_oper_name   = 'adv_z_c2',"                  // &
!                                  "w2p_operator_name       = 'vertical_interp_w2p_sbp42', /", &
!                                  points_type              = "p",                      &
!                                  horizontal_staggering    = "Ah",                     &
!                                  vertical_staggering      = "CharneyPhilips")
! if (parcomm_global%myid==0) then
!     print *, "advection_p_Ah, sbp_d42, adv_z_c2"
!     print "(A,4E25.16)", "Err: ", errs%values
! end if


call deinit_global_parallel_enviroment()

end program test_diffops_3d
