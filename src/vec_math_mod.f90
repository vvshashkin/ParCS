module vec_math_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use parcomm_mod,    only : parcomm_t
use mpi

implicit none

contains

function l2norm(f, mesh, parcomm) result(out)

    type(grid_field_t), intent(in) :: f
    type(mesh_t),       intent(in) :: mesh
    type(parcomm_t),    intent(in) :: parcomm
    real(kind=8)                   :: out

    integer(kind=4) :: t, err
    real(kind=8)    :: out_loc

    out_loc = 0.0_8

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_l2norm_squared_on_tile(f%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

    out = sqrt(out)

end function l2norm

function calc_l2norm_squared_on_tile(f, mesh) result(out)

    type(tile_field_t), intent(in) :: f
    type(tile_mesh_t),  intent(in) :: mesh
    real(kind=8)                   :: out

    integer(kind=4) :: k, j, i

    out = 0.0_8

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                out = out + f%p(i,j,k)**2*mesh%G(i,j)*mesh%hx*mesh%hy
            end do
        end do
    end do

end function calc_l2norm_squared_on_tile


end module vec_math_mod
