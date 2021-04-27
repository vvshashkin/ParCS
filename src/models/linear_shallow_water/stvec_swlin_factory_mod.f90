module stvec_swlin_factory_mod

use stvec_swlin_mod, only : stvec_swlin_t
use partition_mod,   only : partition_t

implicit none

contains

subroutine create_stvec_swlin(new_stvec, halo_width, partition)

    use grid_field_factory_mod, only : create_grid_field

    type(stvec_swlin_t), intent(out) :: new_stvec
    type(partition_t),   intent(in)  :: partition
    integer(kind=4),     intent(in)  :: halo_width

    integer(kind=4) :: t

    call create_grid_field(new_stvec%h, halo_width, 0, partition)
    call create_grid_field(new_stvec%u, halo_width, 0, partition)
    call create_grid_field(new_stvec%v, halo_width, 0, partition)

    new_stvec%ts = partition%ts
    new_stvec%te = partition%te

    do t = partition%ts, partition%te
        new_stvec%h%block(t)%p(:,:,:) = 0._8
        new_stvec%u%block(t)%p(:,:,:) = 0._8
        new_stvec%v%block(t)%p(:,:,:) = 0._8
    end do

    new_stvec%init_and_alloc = .true.

end subroutine create_stvec_swlin

end module stvec_swlin_factory_mod
