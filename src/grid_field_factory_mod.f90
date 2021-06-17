module grid_field_factory_mod

use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t

implicit none

contains

subroutine create_grid_field(grid_field, halo_width_xy, halo_width_z, mesh)

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z

    type(mesh_t),  intent(in)  :: mesh

    integer(kind=4) :: t

    allocate(grid_field%tile(mesh%ts:mesh%te))

    !grid_field%ts = mesh%ts
    !grid_field%te = mesh%te

    do t = mesh%ts, mesh%te
        call grid_field%tile(t)%init(mesh%tile(t)%is, mesh%tile(t)%ie, &
                                     mesh%tile(t)%js, mesh%tile(t)%je, &
                                     mesh%tile(t)%ks, mesh%tile(t)%ke, &
                                     halo_width_xy, halo_width_xy, halo_width_z)
    end do

end subroutine create_grid_field

end module grid_field_factory_mod
