module vec_math_mod

use grid_field_mod, only : grid_field_t, tile_field_t
use mesh_mod,       only : mesh_t, tile_mesh_t
use parcomm_mod,    only : parcomm_t
use mpi

implicit none

contains

function mass(f, mesh, parcomm) result(out)

    type(grid_field_t), intent(in) :: f
    type(mesh_t),       intent(in) :: mesh
    type(parcomm_t),    intent(in) :: parcomm
    real(kind=8)                   :: out

    integer(kind=4) :: t, err
    real(kind=8)    :: out_loc

    out_loc = 0.0_8

    do t = mesh%ts, mesh%te
        out_loc = out_loc + calc_mass_tile(f%tile(t), mesh%tile(t))
    end do

    call mpi_allreduce(out_loc, out, 1, mpi_double, mpi_sum, parcomm%comm_w, err)

end function mass
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

subroutine cube2cart_vec(vx, vy, vz, u, v, mesh)

    type(grid_field_t), intent(in)    :: u, v
    type(grid_field_t), intent(inout) :: vx, vy, vz
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call calc_cube2cart_vec_tile(vx%tile(t), vy%tile(t), vz%tile(t), &
                                     u%tile(t), v%tile(t), mesh%tile(t) )
    end do

end subroutine cube2cart_vec

subroutine calc_cube2cart_vec_tile(vx, vy, vz, u, v, mesh)

    type(tile_field_t), intent(in)    :: u, v
    type(tile_field_t), intent(inout) :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    real(kind=8)    :: v_xyz(3)
    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                v_xyz(1:3) = mesh%b1(1:3,i,j)*u%p(i,j,k)+mesh%b2(1:3,i,j)*v%p(i,j,k)
                vx%p(i,j,k) = v_xyz(1)
                vy%p(i,j,k) = v_xyz(2)
                vz%p(i,j,k) = v_xyz(3)
            end do
        end do
    end do
end subroutine calc_cube2cart_vec_tile

subroutine cart2cube_vec(u, v, vx, vy, vz, mesh)

    type(grid_field_t), intent(inout) :: u, v
    type(grid_field_t), intent(in)    :: vx, vy, vz
    type(mesh_t),       intent(in)    :: mesh

    integer(kind=4) :: t

    do t = mesh%ts, mesh%te
        call calc_cart2cube_vec_tile(u%tile(t), v%tile(t), &
                                     vx%tile(t), vy%tile(t), vz%tile(t), &
                                     mesh%tile(t) )
    end do

end subroutine cart2cube_vec

subroutine calc_cart2cube_vec_tile(u, v, vx, vy, vz, mesh)

    type(tile_field_t), intent(inout) :: u, v
    type(tile_field_t), intent(in)    :: vx, vy, vz
    type(tile_mesh_t),  intent(in)    :: mesh

    integer(kind=4) :: k, j, i

    do k = mesh%ks, mesh%ke
        do j = mesh%js, mesh%je
            do i = mesh%is, mesh%ie
                u%p(i,j,k) = vx%p(i,j,k)*mesh%b1(1,i,j) + &
                             vy%p(i,j,k)*mesh%b1(2,i,j) + &
                             vz%p(i,j,k)*mesh%b1(3,i,j)

                v%p(i,j,k) = vx%p(i,j,k)*mesh%b2(1,i,j) + &
                             vy%p(i,j,k)*mesh%b2(2,i,j) + &
                             vz%p(i,j,k)*mesh%b2(3,i,j)
            end do
        end do
    end do
end subroutine calc_cart2cube_vec_tile

end module vec_math_mod
