program test_diffops_3d

use parcomm_mod,         only : init_global_parallel_enviroment, &
                                deinit_global_parallel_enviroment, parcomm_global
use test_diffops_3d_mod, only : test_scalar_advection_3d
use key_value_mod,       only : key_value_r8_t

implicit none

real(kind=8) :: err
type(key_value_r8_t)  :: errs
integer(kind=4), parameter :: Ns(3) = [32,64,128]

call init_global_parallel_enviroment()

errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
                                 advection_oper_name      = "advection_w_staggered", &
                                 hor_advection_oper_name  = "up4",                   &
                                 vert_advection_oper_name = "adv_z_c2",              &
                                 points_type              = "w",                     &
                                 horizontal_staggering    = "C",                     &
                                 vertical_staggering      = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "advection_w_staggered, up4, adv_z_c2"
    print "(A,4E25.16)", "Err: ", errs%values
end if

errs =  test_scalar_advection_3d(Nh = 32, Nz = 10, &
                                 advection_oper_name      = "advection_p_staggered", &
                                 hor_advection_oper_name  = "up4",                   &
                                 vert_advection_oper_name = "adv_z_c2",              &
                                 points_type              = "p",                     &
                                 horizontal_staggering    = "C",                     &
                                 vertical_staggering      = "CharneyPhilips")
if (parcomm_global%myid==0) then
    print *, "advection_p_staggered, up4, adv_z_c2"
    print "(A,4E25.16)", "Err: ", errs%values
end if
call deinit_global_parallel_enviroment()

end program test_diffops_3d
