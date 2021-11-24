module interpolator_w2uv_factory_mod

use abstract_interpolators3d_mod,  only : interpolator_w2uv_t
use interpolators_w2uv_mod,        only : w2uv_colocated_t, w2uv_hor_colocated_t,&
                                          w2uv_staggered_t
use vertical_operator_factory_mod, only : create_vertical_operator
use parcomm_mod,                   only : parcomm_global


implicit none

contains

subroutine create_w2uv_interpolator(w2uv_interpolator, w2uv_interpolator_name)
    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name

    if(w2uv_interpolator_name == "w2uv_colocated") then
        w2uv_interpolator = w2uv_colocated_t()
    else if(w2uv_interpolator_name(1:18) == "w2uv_hor_colocated") then
        call create_hor_colocated_w2uv(w2uv_interpolator, w2uv_interpolator_name)
    else if(w2uv_interpolator_name(1:14) == "w2uv_staggered") then
        call parcomm_global%abort("not implemented")
    else
        call parcomm_global%abort("create_w2uv_interpolator error, unknown w2uv_interpolator_name: "//&
                                   w2uv_interpolator_name)
    end if
end subroutine create_w2uv_interpolator

subroutine create_hor_colocated_w2uv(w2uv_interpolator, w2uv_interpolator_name)

    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name
    type(w2uv_hor_colocated_t), allocatable :: w2uv_hor_colocated

    allocate(w2uv_hor_colocated)

    select case(w2uv_interpolator_name)
    case("w2uv_hor_colocated_sbp21")
        call create_vertical_operator(w2uv_hor_colocated%w2p_oper,&
                                      "vertical_interp_w2p_sbp21")
    case default
        call parcomm_global%abort("create_hor_colocated_w2uv error, unknown w2uv_interpolator_name: "//&
                                   w2uv_interpolator_name)
    end select

    call move_alloc(w2uv_hor_colocated, w2uv_interpolator)
end subroutine create_hor_colocated_w2uv

end module interpolator_w2uv_factory_mod
