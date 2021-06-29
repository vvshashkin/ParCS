module topology_factory_mod

use topology_mod, only : topology_t

contains

function create_topology(topology_type) result(topology)

    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t

    character(len=*),  intent(in)  :: topology_type
    class(topology_t), allocatable :: topology

    select case(topology_type)
    case("cube")
        allocate(cubed_sphere_topology_t :: topology)
        call topology%init()
    case default
        call avost("unknown topology type " //  topology_type)
    end select
end function create_topology

end module topology_factory_mod
