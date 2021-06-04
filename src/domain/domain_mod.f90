module domain_mod

use topology_mod,  only : topology_t
use partition_mod, only : partition_t
use mesh_mod,      only : mesh_t
use parcomm_mod,   only : parcomm_t

implicit none

type, public :: domain_t

    class(topology_t), allocatable  :: topology
    type(parcomm_t)   :: parcomm
    type(partition_t) :: partition
    type(mesh_t)      :: mesh_u, mesh_v, mesh_p

contains

end type domain_t

contains


end module domain_mod