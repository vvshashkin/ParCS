module default_quadrature_mod

use parcomm_mod,    only : parcomm_t
use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use abstract_quadrature_mod, only : quadrature_t

use mpi

implicit none

type, extends(quadrature_t) :: default_quadrature_t
    contains
        procedure, public :: mass => calc_mass_default
        procedure, public :: dot => calc_dot_default
end type default_quadrature_t

contains

function calc_mass_default(this, f, mesh, parcomm) result(mass)
    class(default_quadrature_t), intent(in) :: this
    type(grid_field_t),          intent(in) :: f
    type(mesh_t),                intent(in) :: mesh
    type(parcomm_t),             intent(in) :: parcomm

    real(kind=8) :: mass

    real(kind=8) :: mass_loc
    integer(kind=4) :: t, err

    mass_loc = 0.0_8

    do t = mesh%ts, mesh%te
        mass_loc = mass_loc + calc_mass_tile(f%tile(t), mesh%tile(t))
    end do
    call mpi_allreduce(mass_loc, mass, 1, mpi_double, mpi_sum, parcomm%comm_w, err)
end

function calc_dot_default(this, f1, f2, mesh, parcomm) result(dot_prod)
    class(default_quadrature_t), intent(in) :: this
    type(grid_field_t),          intent(in) :: f1, f2
    type(mesh_t),                intent(in) :: mesh
    type(parcomm_t),             intent(in) :: parcomm

    real(kind=8) :: dot_prod

    real(kind=8) :: dot_prod_loc
    integer(kind=4) :: t, err

    dot_prod_loc = 0.0_8

    do t = mesh%ts, mesh%te
        dot_prod_loc = dot_prod_loc + calc_dot_tile(f1%tile(t), f2%tile(t), mesh%tile(t))
    end do
    call mpi_allreduce(dot_prod_loc, dot_prod, 1, mpi_double, mpi_sum, parcomm%comm_w, err)
end

function calc_mass_tile(f, mesh) result(out)

    type(tile_field_t), intent(in) :: f
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f%p(i,j,k)*mesh%G(i,j)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_mass_tile

function calc_dot_tile(f1, f2, mesh) result(out)

    type(tile_field_t), intent(in) :: f1, f2
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f1%p(i,j,k)*f2%p(i,j,k)*mesh%G(i,j)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_dot_tile

end module default_quadrature_mod
