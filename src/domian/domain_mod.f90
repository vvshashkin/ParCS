module domain_mod

use partition_mod, only : partition_t
use mesh_mod,      only : mesh_t

implicit none

type, public :: domain_t

    type(partition_t)         :: partition
    type(mesh_t), allocatable :: mesh(:), mesh_u(:), mesh_v(:), mesh_p(:)

contains

end type domain_t

contains


end module domain_mod
