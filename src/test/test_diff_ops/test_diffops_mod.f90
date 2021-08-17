module test_diffops_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
    
implicit none

contains

real(kind=8) function test_div_a2(N) result(err)

    use test_fields_mod,  only : set_vector_test_field, solid_rot=>solid_rotation_field_generator
    use div_factory_mod,  only : create_div_operator
    use abstract_div_mod, only : div_operator_t 

    integer(kind=4), intent(in) :: N
    !locals:
    integer(kind=4), parameter  :: nz = 3
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: u, v, div
    type(domain_t)              :: domain
    class(div_operator_t), allocatable :: div_op

    call create_domain(domain, "cube", 'A', N, nz)
    call create_grid_field(u, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(div, 0, 0, domain%mesh_v)

    call set_vector_test_field(u,v,solid_rot, domain%mesh_u, domain%mesh_v, &
                               0, "contravariant")


    div_op = create_div_operator(domain, "divergence_a2_ecs")
    call div_op%calc_div(div,u,v,domain)
    err = div%maxabs(domain%mesh_p,domain%parcomm)
end function test_div_a2

real(kind=8) function test_grad_a2(N) result(err)

    use test_fields_mod,    only : set_vector_test_field, set_scalar_test_field, &
                                   xyz_f => xyz_scalar_field_generator, &
                                   xyz_grad => xyz_grad_generator
    use grad_factory_mod,   only : create_grad_operator
    use abstract_grad_mod,  only : grad_operator_t 

    integer(kind=4), intent(in) :: N
    !locals:
    integer(kind=4), parameter  :: nz = 3
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: gx, gy, f
    type(grid_field_t)          :: gx_true, gy_true
    type(domain_t)              :: domain
    class(grad_operator_t), allocatable :: grad_op

    call create_domain(domain, "cube", 'A', N, nz)
    call create_grid_field(gx, 0, 0, domain%mesh_u)
    call create_grid_field(gy, 0, 0, domain%mesh_v)
    call create_grid_field(gx_true, 0, 0, domain%mesh_u)
    call create_grid_field(gy_true, 0, 0, domain%mesh_v)
    call create_grid_field(f, ex_halo_width, 0, domain%mesh_v)

    call set_scalar_test_field(f,xyz_f, domain%mesh_p,0)
    call set_vector_test_field(gx_true, gy_true, xyz_grad, domain%mesh_u, domain%mesh_v, 0, "contravariant")

    grad_op = create_grad_operator(domain, "gradient_a2_ecs")
    call grad_op%calc_grad(gx,gy,f,domain)
    call gx%update(-1.0_8,gx_true,domain%mesh_u)
    call gy%update(-1.0_8,gy_true,domain%mesh_v)
    err = gx%maxabs(domain%mesh_p,domain%parcomm)+gy%maxabs(domain%mesh_p,domain%parcomm)
end function test_grad_a2

end module test_diffops_mod