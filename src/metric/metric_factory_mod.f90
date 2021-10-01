module metric_factory_mod

use metric_mod,        only : metric_t
use topology_mod,      only : topology_t
use parcomm_mod,       only : parcomm_global
use config_metric_mod, only : config_metric_t

implicit none

contains

subroutine create_metric_by_config(metric, topology, metric_type, config)
    use ecs_metric_mod,            only : ecs_metric_t
    use ecs_metric_factory_mod,    only : create_ecs_metric

    class(metric_t), allocatable, intent(out) :: metric
    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    type(config_metric_t),        intent(in)  :: config

    select case(metric_type)
    case("ecs")
        call create_ecs_metric(metric, topology, &
         config%scale, config%omega, config%rotation_matrix, config%rotation_axis)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end subroutine create_metric_by_config

subroutine create_metric(metric, topology, metric_type)
    use ecs_metric_mod,            only : ecs_metric_t
    use ecs_metric_factory_mod,    only : create_ecs_metric

    class(topology_t),            intent(in)  :: topology
    character(len=*),             intent(in)  :: metric_type
    class(metric_t), allocatable, intent(out) :: metric

    select case(metric_type)
    case("ecs")
        call create_ecs_metric(metric, topology)
    case default
        call parcomm_global%abort("Unknown metric_type var " // metric_type // &
                                                    " in metric_factory_mod")
    end select
end subroutine create_metric

end module metric_factory_mod
