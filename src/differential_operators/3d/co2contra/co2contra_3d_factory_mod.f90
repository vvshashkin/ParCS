module co2contra_3d_factory_mod

use domain_mod,                 only : domain_t
use abstract_co2contra_3d_mod,  only : co2contra_3d_operator_t
use co2contra_3d_colocated_mod, only : co2contra_3d_colocated_t
use parcomm_mod,                only : parcomm_global

implicit none

contains

subroutine create_co2contra_3d_operator(co2contra_3d_op, domain, co2contra_3d_name)
    class(co2contra_3d_operator_t), allocatable , intent(out) :: co2contra_3d_op
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in) :: co2contra_3d_name

    if (co2contra_3d_name == "co2contra_3d_colocated") then
        co2contra_3d_op = co2contra_3d_colocated_t()
    else if (co2contra_3d_name(1:19) == "co2contra_3d_Cgrid_") then
        call create_co2contra_3d_Cgrid(co2contra_3d_op, domain, co2contra_3d_name(20:))
    else if (co2contra_3d_name(1:25) == "co2contra_3d_h_colocated_") then
        call create_co2contra_3d_h_colocated(co2contra_3d_op, domain, co2contra_3d_name(26:))
    else
        call parcomm_global%abort("unknown co2contra_3d operator: "//co2contra_3d_name)
    end if

end subroutine create_co2contra_3d_operator

subroutine create_co2contra_3d_Cgrid(co2contra_3d_op, domain, co2contra_3d_name)

    use co2contra_3d_Cgrid_mod,        only : co2contra_3d_Cgrid_t
    use vertical_operator_factory_mod, only : create_vertical_operator
    use interpolator_h2v_factory_mod,  only : create_h2v_interpolator
    use interpolator_v2h_factory_mod,  only : create_v2h_interpolator
    use grid_field_factory_mod,        only : create_grid_field


    class(co2contra_3d_operator_t), allocatable , intent(out) :: co2contra_3d_op
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in) :: co2contra_3d_name

    integer(kind=4) :: halo_width

    type(co2contra_3d_Cgrid_t), allocatable :: co2contra

    allocate(co2contra)

    if (co2contra_3d_name(1:8) == "h_sbp21_") then
        halo_width = 4
        call create_h2v_interpolator(co2contra%interp_h2v, "W21_stagered_interp_c2i", domain)
        call create_v2h_interpolator(co2contra%interp_v2h, "W21_stagered_interp_i2c", domain)
    else if (co2contra_3d_name(1:8) == "h_sbp42_") then
        halo_width = 4
        call create_h2v_interpolator(co2contra%interp_h2v, "W42_stagered_interp_c2i", domain)
        call create_v2h_interpolator(co2contra%interp_v2h, "W42_stagered_interp_i2c", domain)
    else
        call parcomm_global%abort("unknown co2contra_3d_Cgrid operator: "//co2contra_3d_name)
    end if
    if (co2contra_3d_name(9:) == "v_sbp21") then
        call create_vertical_operator(co2contra%interp_p2w, "vertical_interp_p2w_sbp21")
        call create_vertical_operator(co2contra%interp_w2p, "vertical_interp_w2p_sbp21")
    else if (co2contra_3d_name(9:) == "v_sbp42") then
        call create_vertical_operator(co2contra%interp_p2w, "vertical_interp_p2w_sbp42")
        call create_vertical_operator(co2contra%interp_w2p, "vertical_interp_w2p_sbp42")
    else if (co2contra_3d_name(9:) == "v_colocated") then
        call create_vertical_operator(co2contra%interp_p2w, "identity_p")
        call create_vertical_operator(co2contra%interp_w2p, "identity_p")
    else
        call parcomm_global%abort("unknown co2contra_3d_Cgrid operator: "//co2contra_3d_name)
    end if

    call create_grid_field(co2contra%up, halo_width, 0, domain%mesh_u)
    call create_grid_field(co2contra%vp, halo_width, 0, domain%mesh_v)
    call create_grid_field(co2contra%wp, 0,          0, domain%mesh_w)

    call move_alloc(co2contra, co2contra_3d_op)

end subroutine create_co2contra_3d_Cgrid

subroutine create_co2contra_3d_h_colocated(co2contra_3d_op, domain, co2contra_3d_name)

    use co2contra_3d_h_colocated_mod,  only : co2contra_3d_h_colocated_t
    use vertical_operator_factory_mod, only : create_vertical_operator
    use grid_field_factory_mod,        only : create_grid_field


    class(co2contra_3d_operator_t), allocatable , intent(out) :: co2contra_3d_op
    type(domain_t),   intent(in) :: domain
    character(len=*), intent(in) :: co2contra_3d_name

    integer(kind=4) :: halo_width

    type(co2contra_3d_h_colocated_t), allocatable :: co2contra

    allocate(co2contra)

    if (co2contra_3d_name == "v_sbp21") then
        call create_vertical_operator(co2contra%interp_p2w, "vertical_interp_p2w_sbp21")
        call create_vertical_operator(co2contra%interp_w2p, "vertical_interp_w2p_sbp21")
    else if (co2contra_3d_name == "v_sbp42") then
        call create_vertical_operator(co2contra%interp_p2w, "vertical_interp_p2w_sbp42")
        call create_vertical_operator(co2contra%interp_w2p, "vertical_interp_w2p_sbp42")
    !Unecessary option
    else if (co2contra_3d_name == "v_colocated") then
        call create_vertical_operator(co2contra%interp_p2w, "identity_p")
        call create_vertical_operator(co2contra%interp_w2p, "identity_p")
    else
        call parcomm_global%abort("unknown co2contra_3d_h_colocated operator: "//co2contra_3d_name)
    end if

    call create_grid_field(co2contra%wp, 0, 0, domain%mesh_w)

    call move_alloc(co2contra, co2contra_3d_op)

end subroutine create_co2contra_3d_h_colocated

end module co2contra_3d_factory_mod
