module test_metric_class_mod

implicit none

real(kind=8), parameter :: test_tolerace = 3e-15_8

contains

subroutine test_metric_class(topology_type,metric_type)
    use metric_mod,           only : metric_t
    use metric_factory_mod,   only : create_metric
    use topology_mod,         only : topology_t
    use topology_factory_mod, only : init_topology

    character(len=*), intent(in) :: topology_type, metric_type

    class(topology_t), allocatable :: topology
    class(metric_t),   allocatable :: metric
    real(kind=8) a1(3), a2(3), b1(3), b2(3), r(3)
    real(kind=8) Q(3), QI(3), G
    real(kind=8), parameter :: step = 0.1
    integer(kind=4), parameter :: N = 1/step
    integer(kind=4) npanels, panel_ind, i, j
    real(kind=8) x, y, x0, y0, dx, dy
    logical :: is_correct = .true.

    topology = init_topology(topology_type)
    call create_metric(topology, metric_type, metric)

    npanels = topology%npanels
    x0 = metric%x0
    y0 = metric%y0
    dx = metric%x1-metric%x0
    dx = metric%y1-metric%y0

    do panel_ind = 1, npanels
        do j=0,N
            do i = 0,N
                x = x0+i*step
                y = y0+j*step
                a1 = metric%a1(panel_ind,x,y)
                a2 = metric%a2(panel_ind,x,y)
                b1 = metric%b1(panel_ind,x,y)
                b2 = metric%b2(panel_ind,x,y)
                Q  = metric%Q(panel_ind,x,y)
                QI = metric%QI(panel_ind,x,y)
                G  = metric%G(panel_ind,x,y)
                call check("a1 is not normal to b2",sum(a1(1:3)*b2(1:3)), 0.0_8, test_tolerace, is_correct)
                call check("a2 is not normal to b1",sum(a2(1:3)*b1(1:3)), 0.0_8, test_tolerace, is_correct)
                call check("Q(1) /= a1*a1",sum(a1(1:3)*a1(1:3)), Q(1), test_tolerace, is_correct)
                call check("Q(2) /= a1*a2",sum(a1(1:3)*a2(1:3)), Q(2), test_tolerace, is_correct)
                call check("Q(3) /= a2*a2",sum(a2(1:3)*a2(1:3)), Q(3), test_tolerace, is_correct)

                call check("QI(1) /= b1*b1",sum(b1(1:3)*b1(1:3)), QI(1), test_tolerace, is_correct)
                call check("QI(2) /= b1*b2",sum(b1(1:3)*b2(1:3)), QI(2), test_tolerace, is_correct)
                call check("QI(3) /= b2*b2",sum(b2(1:3)*b2(1:3)), QI(3), test_tolerace, is_correct)

                call check("G**2 /= Q(1)*Q(3)-Q(2)**2",G**2, Q(1)*Q(3)-Q(2)**2, test_tolerace, is_correct)

                if(metric_type == "ecs") then
                    r = metric%point_r(panel_ind,x,y)
                    call check("a1 is not normal to r",sum(a1(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                    call check("a2 is not normal to r",sum(a2(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                    call check("b1 is not normal to r",sum(b1(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                    call check("b2 is not normal to r",sum(b2(1:3)*r(1:3)), 0.0_8, test_tolerace, is_correct)
                end if

                if(.not. is_correct) then
                    print *, "metric_class test failed"
                    return
                end if
            end do
        end do
    end do
    print *, "metric_class test passed"
end subroutine test_metric_class

subroutine check(message,a,b,tol, is_correct)
    !check if abs(a-b) < tol and prints message otherwise
    character(len=*), intent(in)    :: message
    real(kind=8),     intent(in)    :: a, b
    real(kind=8),     intent(in)    :: tol
    logical,          intent(inout) :: is_correct

    logical is_correct_local

    is_correct_local = abs(a-b) < tol
    is_correct = is_correct .and. is_correct_local

    if(.not. is_correct_local) then
        print *, message,", error: ", a-b
    end if

end

end module test_metric_class_mod
