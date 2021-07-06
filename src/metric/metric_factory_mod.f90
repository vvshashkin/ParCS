module metric_factory_mod

use metric_mod,   only : metric_t
use topology_mod, only : topology_t
use parcomm_mod,  only : parcomm_global

implicit none

contains

subroutine create_metric(metric,topology,metric_type)
    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
    use ecs_metric_mod,            only : ecs_metric_t
    use ecs_metric_factory_mod,    only : create_ecs_metric

    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    class(metric_t), allocatable, intent(out) :: metric

    select case(metric_type)
    case("ecs")
        select type(topology)
        class is (cubed_sphere_topology_t)
            call create_ecs_metric(metric,topology)
        class default
            call parcomm_global%abort( &
                 "Wrong topology class in metric_factory_mod. ecs case")
        end select
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end

end
