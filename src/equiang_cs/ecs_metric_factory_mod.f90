module ecs_metric_factory_mod

implicit none

contains

subroutine create_ecs_metric(topology, metric, sphere_r, rotation_matrix)
    use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
    use metric_mod,                only : metric_t
    use ecs_metric_mod,            only : ecs_metric_t
    use const_mod,                 only : radz, pi

    class(cubed_sphere_topology_t), intent(in)    :: topology
    class(metric_t), allocatable,   intent(inout) :: metric
    real(kind=8), optional,         intent(in)    :: sphere_r
    real(kind=8), optional,         intent(in)    :: rotation_matrix(3,3)
    !local
    real(kind=8) rotation_matrix_local(3,3)

    allocate(ecs_metric_t :: metric)
    metric%x0 =-0.25*pi
    metric%y0 =-0.25*pi
    metric%x1 = 0.25*pi
    metric%y1 = 0.25*pi

    select type(metric)
    class is (ecs_metric_t)
        metric%topology = topology

        if(present(sphere_r)) then
            metric%a = sphere_r
        else
            metric%a = radz
        end if

        if(present(rotation_matrix)) then
            metric%rotation_matrix = rotation_matrix
        else
            metric%rotation_matrix(1:3,1:3) = reshape((/1.0_8, 0.0_8, 0.0_8,  &
                                                        0.0_8, 1.0_8, 0.0_8,  &
                                                        0.0_8, 0.0_8, 1.0_8/),&
                                                        (/3,3/))
        end if
    class default
        call avost("strange mistake in create_ecs_metric")
    end select
end

end module ecs_metric_factory_mod
