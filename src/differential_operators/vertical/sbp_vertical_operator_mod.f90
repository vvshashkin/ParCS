module sbp_vertical_operator_mod

use abstract_vertical_operator_mod, only : vertical_operator_t
use domain_mod,                     only : domain_t
use grid_field_mod,                 only : grid_field_t
use sbp_operator_mod,               only : sbp_operator_t
use parcomm_mod,                    only : parcomm_global

implicit none

type, extends(vertical_operator_t) :: sbp_vertical_op_t
    class(sbp_operator_t), allocatable :: sbp_op
    character(len=:),      allocatable :: tagret_grid
    contains
    procedure :: apply
end type sbp_vertical_op_t

contains

subroutine apply(this, f_out, f_in, domain)
    class(sbp_vertical_op_t), intent(in)    :: this
    type(grid_field_t),       intent(inout) :: f_out
    type(grid_field_t),       intent(in)    :: f_in
    type(domain_t),           intent(in)    :: domain

    real(kind=8) :: scale, hz

    !call f_out%assign(0.0_8,domain%mesh_n)
    select case(this%tagret_grid)
    case("n")
        call this%sbp_op%apply_z(f_out,f_in,domain%mesh_n)
        scale = domain%mesh_n%vertical_scale
        hz    = domain%mesh_n%tile(domain%mesh_n%ts)%hz
        call f_out%assign(1.0_8/(scale*hz), f_out, domain%mesh_n)
    case("p")
        call this%sbp_op%apply_z(f_out,f_in,domain%mesh_p)
        scale = domain%mesh_p%vertical_scale
        hz    = domain%mesh_p%tile(domain%mesh_p%ts)%hz
        call f_out%assign(1.0_8/(scale*hz), f_out, domain%mesh_p)
    case default
        call parcomm_global%abort("vertical sbp operator apply error "// &
                                  "unsupported target grid: "// this%tagret_grid)
    end select

end subroutine apply

end module sbp_vertical_operator_mod
