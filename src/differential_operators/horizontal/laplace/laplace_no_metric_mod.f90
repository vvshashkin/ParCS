module laplace_no_metric_mod

use abstract_laplace_mod,   only : laplace_operator_t
use grid_field_mod,         only : grid_field_t
use domain_mod,             only : domain_t
use sbp_operator_mod,       only : sbp_operator_t
use exchange_abstract_mod,  only : exchange_t
use halo_mod,               only : halo_t

implicit none

type, extends(laplace_operator_t) :: laplace_no_metric_t
    type(grid_field_t)             :: d2f_x, d2f_y
    type(sbp_operator_t)           :: sbp_d2_op
    class(halo_t), allocatable     :: edge_sync
    class(exchange_t), allocatable :: exchange
contains
    procedure :: calc_laplace
end type laplace_no_metric_t

contains

subroutine calc_laplace(this, f1, f, domain)

    use vec_math_mod,   only : divide_by_J_self

    class(laplace_no_metric_t), intent(inout) :: this
    type(domain_t),             intent(in)    :: domain
    type(grid_field_t),         intent(inout) :: f
    !output:
    type(grid_field_t),         intent(inout) :: f1

    integer(kind=4) :: t

    call this%exchange%do(f, domain%parcomm)

    do t = domain%partition%ts, domain%partition%te

        call calc_laplace_tile(f1%tile(t), this%d2f_x%tile(t), this%d2f_y%tile(t), f%tile(t), &
                               this%sbp_d2_op, domain%mesh_xy%tile(t), domain%mesh_xy%scale)

    end do

    ! call divide_by_J_self(f1, domain%mesh_xy)

    call this%edge_sync%get_halo_scalar(f1,domain,0)

end subroutine calc_laplace

subroutine calc_laplace_tile(fout, d2f_x, d2f_y, f, sbp_d2_op, mesh_xy, scale)

    use grid_field_mod, only : grid_field_t, tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use tile_mod,       only : tile_t

    !result
    type(tile_field_t), intent(inout) :: fout, d2f_x, d2f_y
    !input
    type(tile_field_t),   intent(in) :: f
    type(sbp_operator_t), intent(in) :: sbp_d2_op
    type(tile_mesh_t),    intent(in) :: mesh_xy
    real(kind=8),         intent(in) :: scale

    type(tile_t) :: work

    work = tile_t(is = mesh_xy%is, ie = mesh_xy%ie, &
                  js = mesh_xy%js, je = mesh_xy%je, &
                  ks = mesh_xy%ks, ke = mesh_xy%ke)

    call sbp_d2_op%apply      (d2f_x, work, mesh_xy%nx, "x", f)
    call sbp_d2_op%add_penalty(d2f_x, work, mesh_xy%nx, "x", "at_interface_new_temp", f)

    call sbp_d2_op%apply      (d2f_y, work, mesh_xy%ny, "y", f)
    call sbp_d2_op%add_penalty(d2f_y, work, mesh_xy%ny, "y", "at_interface_new_temp", f)
    call fout%assign((mesh_xy%hx*scale)**(-2), d2f_x, &
                     (mesh_xy%hy*scale)**(-2), d2f_y, mesh_xy)

end subroutine calc_laplace_tile

end module laplace_no_metric_mod
