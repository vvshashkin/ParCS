module grid_field_factory_mod

use grid_field_mod, only : grid_field_t
use partition_mod,  only : partition_t

implicit none

contains

subroutine create_grid_field(grid_field, halo_width_xy, halo_width_z, partition)

    type(grid_field_t), intent(out) :: grid_field
    integer(kind=4),    intent(in)  :: halo_width_xy, halo_width_z

    type(partition_t),  intent(in)  :: partition

    integer(kind=4) :: t

    allocate(grid_field%block(partition%ts:partition%te))

    do t = partition%ts, partition%te
        call grid_field%block(t)%init(partition%tile(t)%panel_number, &
                        partition%tile(t)%is, partition%tile(t)%ie  , &
                        partition%tile(t)%js, partition%tile(t)%je  , &
                        partition%tile(t)%ks, partition%tile(t)%ke  , &
                        halo_width_xy, halo_width_xy, halo_width_z)
    end do


end subroutine create_grid_field

end module grid_field_factory_mod
