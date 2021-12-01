module interpolator_uv2w_factory_mod

use abstract_interpolators3d_mod,  only : interpolator_uv2w_t
use interpolators_uv2w_mod,        only : uv2w_staggered_t
use vertical_operator_factory_mod, only : create_vertical_operator
use domain_mod,                    only : domain_t
use parcomm_mod,                   only : parcomm_global

implicit none

contains

subroutine create_uv2w_interpolator(uv2w_interpolator, uv2w_interpolator_name, domain)
    class(interpolator_uv2w_t), allocatable, intent(out) :: uv2w_interpolator
    character(len=*), intent(in) :: uv2w_interpolator_name
    type(domain_t),   intent(in) :: domain

    if(uv2w_interpolator_name(1:15) == "uv2w_staggered_") then
        call create_staggered_uv2w(uv2w_interpolator, uv2w_interpolator_name, domain)
    else
        call parcomm_global%abort("create_uv2w_interpolator error, unknown uv2w_interpolator_name: "//&
                                   uv2w_interpolator_name)
    end if
end subroutine create_uv2w_interpolator

subroutine create_staggered_uv2w(uv2w_interpolator, uv2w_interpolator_name, domain)

    use interpolator_v2h_factory_mod, only : create_v2h_interpolator
    use grid_field_factory_mod,       only : create_grid_field

    class(interpolator_uv2w_t), allocatable, intent(out) :: uv2w_interpolator
    character(len=*), intent(in) :: uv2w_interpolator_name
    type(domain_t),   intent(in) :: domain

    type(uv2w_staggered_t), allocatable :: uv2w_staggered

    integer(kind=4) :: halo_width
    integer(kind=4) :: l1 = len("uv2w_staggered_") !operator name prefix name
    integer(kind=4) :: l2

    allocate(uv2w_staggered)

    !sanity check
    if(uv2w_interpolator_name(1:l1) /= "uv2w_staggered_") call local_throw_error()

    !select horizontal part
    if(uv2w_interpolator_name(l1+1:l1+8) == "C_sbp21_") then
        call create_v2h_interpolator(uv2w_staggered%v2p_oper,  &
                                     "W21_stagered_interp_i2c", domain)
        halo_width = 1
        l2 = l1+8 !length of "uv2w_staggered_C_sbp21_"
    else if(uv2w_interpolator_name(l1+1:l1+8) == "C_sbp42_") then
        call create_v2h_interpolator(uv2w_staggered%v2p_oper,  &
                                     "W42_stagered_interp_i2c", domain)
        halo_width = 3
        l2 = l1+8
    else
        call local_throw_error()
    end if

    call create_grid_field(uv2w_staggered%up,0,0,domain%mesh_p)
    call create_grid_field(uv2w_staggered%vp,0,0,domain%mesh_p)

    !select vertical part
    if(uv2w_interpolator_name(l2+1:l2+7) == "v_sbp21") then
        call create_vertical_operator(uv2w_staggered%p2w_oper,&
                                        "vertical_interp_p2w_sbp21")
    else if(uv2w_interpolator_name(l2+1:l2+7) == "v_sbp42") then
        call create_vertical_operator(uv2w_staggered%p2w_oper,&
                                        "vertical_interp_p2w_sbp42")
    else
        call local_throw_error()
    end if

    call move_alloc(uv2w_staggered, uv2w_interpolator)

    contains
    subroutine local_throw_error()
        call parcomm_global%abort("create_staggered_uv2w, unsupported uv2w_interpolator_name "//&
                                  uv2w_interpolator_name)
    end subroutine local_throw_error
end subroutine create_staggered_uv2w

end module interpolator_uv2w_factory_mod
