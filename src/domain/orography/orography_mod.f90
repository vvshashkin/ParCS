module orography_mod

use grid_field_mod, only : grid_field_t

implicit none

type orography_t
    type(orography_1mesh_t), allocatable :: orography_o, orography_x, &
                                            orography_y, orography_xy
end type orography_t

type orography_1mesh_t
    type(grid_field_t) :: h, dh_alpha, dh_beta
end type orography_1mesh_t

end module orography_mod
