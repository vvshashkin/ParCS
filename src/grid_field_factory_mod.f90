module grid_field_factory_mod

use grid_field_mod, only : grid_field_t

implicit none

contains

subroutine create_grid_field(grid_field, halo_width_xy, halo_width_z, mesh)

    use mesh_mod,       only : mesh_t

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z

    type(mesh_t),  intent(in)  :: mesh

    integer(kind=4) :: t

    allocate(grid_field%tile(mesh%ts:mesh%te))

    grid_field%ts = mesh%ts
    grid_field%te = mesh%te

    do t = mesh%ts, mesh%te
        call grid_field%tile(t)%init(mesh%tile(t)%is-halo_width_xy, mesh%tile(t)%ie+halo_width_xy, &
                                     mesh%tile(t)%js-halo_width_xy, mesh%tile(t)%je+halo_width_xy, &
                                     mesh%tile(t)%ks-halo_width_z,  mesh%tile(t)%ke+halo_width_z)
    end do

end subroutine create_grid_field

subroutine create_grid_field_global(grid_field, halo_width_xy, halo_width_z, tile)

    use tile_mod, only : tile_t

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z
    type(tile_t),       intent(in)  :: tile(:)

    integer(kind=4) :: t, ts, te

    ts = lbound(tile, 1)
    te = ubound(tile, 1)

    allocate(grid_field%tile(ts:te))

    grid_field%ts = ts
    grid_field%te = te

    do t = ts, te
        call grid_field%tile(t)%init(tile(t)%is-halo_width_xy, tile(t)%ie+halo_width_xy, &
                                     tile(t)%js-halo_width_xy, tile(t)%je+halo_width_xy, &
                                     tile(t)%ks-halo_width_z,  tile(t)%ke+halo_width_z)
    end do

end subroutine create_grid_field_global

end module grid_field_factory_mod
