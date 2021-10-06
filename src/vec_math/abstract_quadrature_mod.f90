module abstract_quadrature_mod

use parcomm_mod,    only : parcomm_global, parcomm_t
use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t

implicit none

type, abstract :: quadrature_t
    contains
        procedure, public :: mass => calc_mesh_mass_not_implemented
        procedure, public :: dot => calc_mesh_dot_not_implemented
end type quadrature_t

contains

function calc_mesh_mass_not_implemented(this, f, mesh, parcomm) result(mass)
    class(quadrature_t), intent(in) :: this
    type(grid_field_t),  intent(in) :: f
    type(mesh_t),        intent(in) :: mesh
    type(parcomm_t),     intent(in) :: parcomm

    real(kind=8) :: mass

    mass = 0.0_8

    call parcomm_global%abort("tried to use not implemented function mass of quadrature, subtype of mesh_quadrature_t class")
end

function calc_mesh_dot_not_implemented(this, f1, f2, mesh, parcomm) result(dot_prod)
    class(quadrature_t), intent(in) :: this
    type(grid_field_t),  intent(in) :: f1, f2
    type(mesh_t),        intent(in) :: mesh
    type(parcomm_t),     intent(in) :: parcomm

    real(kind=8) :: dot_prod

    dot_prod = 0.0_8

    call parcomm_global%abort("tried to use not implemented function mass of quadrature, subtype of mesh_quadrature_t class")
end

end module abstract_quadrature_mod
