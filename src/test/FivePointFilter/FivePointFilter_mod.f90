module FivePointFilter_mod

use grid_function_mod, only : grid_function_t
use tile_mod,          only : tile_t

implicit none


contains

subroutine calc_fivepointfilter(f1, f2, tile)

    type(grid_function_t), intent(inout) :: f1, f2
    type(tile_t),          intent(in)    :: tile

    integer(kind=4) :: i, j, k

    do k = tile%ks, tile%ke
        do j = tile%js, tile%je
            do i = tile%is, tile%ie
                f2.p(i,j,k) = 0.125_8 * (f1.p(i-1,j,k) + f1.p(i+1,j,k) + f1.p(i,j+1,k) + f1.p(i,j-1,k) + 4.0_8*f1.p(i,j,k))
                print*, ''
                print*,  i, j, k, 0.125_8 * (f1.p(i-1,j,k) + f1.p(i+1,j,k) + f1.p(i,j+1,k) + f1.p(i,j-1,k) + 4.0_8*f1.p(i,j,k))
                print*, ''
            end do
        end do
    end do

end subroutine calc_fivepointfilter

end module FivePointFilter_mod
