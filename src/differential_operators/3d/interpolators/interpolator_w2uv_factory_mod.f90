module interpolator_w2uv_factory_mod

use abstract_interpolators3d_mod,  only : interpolator_w2uv_t
use interpolators_w2uv_mod,        only : w2uv_colocated_t, w2uv_hor_colocated_t,&
                                          w2uv_staggered_t
use vertical_operator_factory_mod, only : create_vertical_operator
use domain_mod,                    only : domain_t
use parcomm_mod,                   only : parcomm_global


implicit none

contains

subroutine create_w2uv_interpolator(w2uv_interpolator, w2uv_interpolator_name, domain)
    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name
    type(domain_t),   intent(in) :: domain

    if(w2uv_interpolator_name == "w2uv_colocated") then
        w2uv_interpolator = w2uv_colocated_t()
    else if(w2uv_interpolator_name(1:18) == "w2uv_hor_colocated") then
        call create_hor_colocated_w2uv(w2uv_interpolator, w2uv_interpolator_name)
    else if(w2uv_interpolator_name(1:15) == "w2uv_staggered_") then
        call create_staggered_w2uv(w2uv_interpolator, w2uv_interpolator_name, domain)
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
    case("w2uv_hor_colocated_sbp42")
        call create_vertical_operator(w2uv_hor_colocated%w2p_oper,&
                                      "vertical_interp_w2p_sbp42")
    case default
        call parcomm_global%abort("create_hor_colocated_w2uv error, unknown w2uv_interpolator_name: "//&
                                   w2uv_interpolator_name)
    end select

    call move_alloc(w2uv_hor_colocated, w2uv_interpolator)
end subroutine create_hor_colocated_w2uv

subroutine create_staggered_w2uv(w2uv_interpolator, w2uv_interpolator_name, domain)

    use interpolator_h2v_factory_mod, only : create_h2v_interpolator
    use grid_field_factory_mod,       only : create_grid_field

    class(interpolator_w2uv_t), allocatable, intent(out) :: w2uv_interpolator
    character(len=*), intent(in) :: w2uv_interpolator_name
    type(domain_t),   intent(in) :: domain

    type(w2uv_staggered_t), allocatable :: w2uv_staggered

    integer(kind=4) :: halo_width
    integer(kind=4) :: l1 = len("w2uv_staggered_") !operator name prefix name
    integer(kind=4) :: l2

    allocate(w2uv_staggered)

    !sanity check
    if(w2uv_interpolator_name(1:l1) /= "w2uv_staggered_") call local_throw_error()

    !select horizontal part
    if(w2uv_interpolator_name(l1+1:l1+8) == "C_sbp21_") then
        call create_h2v_interpolator(w2uv_staggered%h2v_oper,  &
                                     "W21_stagered_interp_c2i", domain)
        halo_width = 1
        l2 = l1+8 !length of "w2uv_staggered_C_sbp21_"
    else if(w2uv_interpolator_name(l1+1:l1+8) == "C_sbp42_") then
        call create_h2v_interpolator(w2uv_staggered%h2v_oper,  &
                                     "W42_stagered_interp_c2i", domain)
        halo_width = 3
        l2 = l1+8
    else
        call local_throw_error()
    end if

    call create_grid_field(w2uv_staggered%wp,halo_width+1,0,domain%mesh_p)

    !select vertical part
    if(w2uv_interpolator_name(l2+1:l2+7) == "v_sbp21") then
        call create_vertical_operator(w2uv_staggered%w2p_oper,&
                                        "vertical_interp_w2p_sbp21")
    else if(w2uv_interpolator_name(l2+1:l2+7) == "v_sbp42") then
        call create_vertical_operator(w2uv_staggered%w2p_oper,&
                                        "vertical_interp_w2p_sbp42")
    else
        call local_throw_error()
    end if

    call move_alloc(w2uv_staggered, w2uv_interpolator)

    contains
    subroutine local_throw_error()
        call parcomm_global%abort("create_staggered_w2uv, unsupported w2uv_interpolator_name "//&
                                  w2uv_interpolator_name)
    end subroutine local_throw_error
end subroutine create_staggered_w2uv

end module interpolator_w2uv_factory_mod
