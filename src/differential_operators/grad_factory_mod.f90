module grad_factory_mod

use abstract_grad_mod, only : grad_operator_t
use domain_mod,        only : domain_t
use parcomm_mod,       only : parcomm_global

implicit none

contains

function create_grad_operator(domain, grad_operator_name) result(grad)
    type(domain_t),    intent(in)  :: domain
    character(len=*),  intent(in)  :: grad_operator_name

    class(grad_operator_t), allocatable :: grad

    if(grad_operator_name == 'gradient_c2_ecs') then
        grad = create_grad_contra_c2_ecs_operator(domain)
    else if(grad_operator_name == 'gradient_c2_cons') then
        grad = create_grad_contra_c2_cons_operator(domain)
    else if(grad_operator_name == 'gradient_c_sbp21') then
        grad = create_grad_contra_c_sbp21_operator(domain)
    else if(grad_operator_name == 'gradient_a2_ecs' .or. &
            grad_operator_name == 'gradient_a2_cons') then
        grad = create_grad_contra_a2_operator(domain,grad_operator_name)
    else if(grad_operator_name == 'gradient_ah2_ecs') then
        grad = create_grad_contra_ah2_operator(domain, grad_operator_name)
    else if(grad_operator_name == 'gradient_ah42_sbp_ecs' .or. &
            grad_operator_name == 'gradient_ah43_sbp_ecs') then
        grad = create_grad_contra_ah_sbp_operator(domain, grad_operator_name)
    else
        call parcomm_global%abort("unknown gradient operator: "//grad_operator_name)
    end if
end

function create_grad_contra_c2_ecs_operator(domain) result(grad)

    use grad_contra_c2_ecs_mod, only : grad_contra_c2_ecs_t
    use halo_factory_mod,       only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(grad_contra_c2_ecs_t)    :: grad

    integer(kind=4), parameter :: ecs_halo_width=2

    grad = grad_contra_c2_ecs_t()
    call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")

end function create_grad_contra_c2_ecs_operator

function create_grad_contra_c_sbp21_operator(domain) result(grad)

    use grad_contra_c2_ecs_mod,   only : grad_contra_c_sbp21_t
    use exchange_factory_mod,     only : create_symm_halo_exchange_A

    type(domain_t),   intent(in)      :: domain
    type(grad_contra_c_sbp21_t)       :: grad

    integer(kind=4), parameter :: halo_width=2

    grad%exch_halo = create_symm_halo_exchange_A( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')
end function create_grad_contra_c_sbp21_operator

function create_grad_contra_c2_cons_operator(domain) result(grad)

    use grad_contra_c2_ecs_mod, only : grad_contra_c2_cons_t
    use halo_factory_mod,   only : create_halo_procedure, create_vector_halo_procedure

    type(domain_t),   intent(in)  :: domain
    type(grad_contra_c2_cons_t)   :: grad

    integer(kind=4), parameter :: halo_width=1

    grad = grad_contra_c2_cons_t()
    call create_halo_procedure(grad%halo_procedure,domain,halo_width,"A_default")
    call create_vector_halo_procedure(grad%sync_procedure,domain,0,"C_vec_default")

end function create_grad_contra_c2_cons_operator

function create_grad_contra_a2_operator(domain, grad_operator_name) result(grad)

    use grad_contra_a2_mod, only : grad_contra_a2_t
    use halo_factory_mod,   only : create_halo_procedure

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_contra_a2_t) :: grad

    integer(kind=4), parameter :: ecs_halo_width=2, default_halo_width=1

    grad = grad_contra_a2_t()
    if(grad_operator_name=="gradient_a2_ecs") then
        call create_halo_procedure(grad%halo_procedure,domain,ecs_halo_width,"ECS_O")
    else if(grad_operator_name=="gradient_a2_cons") then
        call create_halo_procedure(grad%halo_procedure,domain,default_halo_width,"A_default")
    else
        call parcomm_global%abort("grad_factory_mod, create_grad_contra_a2_operator "//&
                                  "unknown gradient_a2 subtype: "// grad_operator_name)
    end if

end function create_grad_contra_a2_operator

function create_grad_contra_ah2_operator(domain, grad_operator_name) result(grad)

    use grad_contra_ah2_mod,  only : grad_contra_ah2_t
    use exchange_factory_mod, only : create_symm_halo_exchange_Ah
    use ecs_metric_mod,       only : ecs_b1_proto, ecs_b2_proto, ecs_a1_proto, ecs_a2_proto
    use const_mod,            only : pi

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_contra_ah2_t)       :: grad

    integer(kind=4), parameter :: halo_width=2

    integer(kind=4) :: t, is, ie, js, je, i, j
    real(kind=8)    :: alpha, beta, a(3), b(3)

    grad%exch_halo = create_symm_halo_exchange_Ah( &
                    domain%partition, domain%parcomm, domain%topology,  halo_width, 'full')
    grad%subtype = grad_operator_name

    allocate(grad%q(domain%mesh_xy%ts : domain%mesh_xy%te))
    do t=domain%mesh_xy%ts, domain%mesh_xy%te
        is = domain%mesh_xy%tile(t)%is
        ie = domain%mesh_xy%tile(t)%ie
        js = domain%mesh_xy%tile(t)%js
        je = domain%mesh_xy%tile(t)%je

        if(js == 1) then
            allocate(grad%q(t)%qb(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a2_proto(alpha,-0.25_8*pi)
                b(1:3) = ecs_b1_proto(alpha, 0.25_8*pi)
                b(1:3) = [b(1),-b(3),b(2)]
                grad%q(t)%qb(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(je == domain%mesh_xy%tile(t)%ny+1) then
            allocate(grad%q(t)%qt(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a2_proto(alpha, 0.25_8*pi)
                b(1:3) = ecs_b1_proto(alpha,-0.25_8*pi)
                b(1:3) = [b(1),b(3),-b(2)]
                grad%q(t)%qt(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(is == 1) then
            allocate(grad%q(t)%ql(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a1_proto(-0.25_8*pi, beta)
                b(1:3) = ecs_b2_proto( 0.25_8*pi, beta)
                b(1:3) = [-b(3),b(2),b(1)]
                grad%q(t)%ql(j) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(ie == domain%mesh_xy%tile(t)%nx+1) then
            allocate(grad%q(t)%qr(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                a(1:3) = ecs_a1_proto( 0.25_8*pi, beta)
                b(1:3) = ecs_b2_proto(-0.25_8*pi, beta)
                b(1:3) = [b(3),b(2),-b(1)]
                grad%q(t)%qr(j) = sum(a(1:3)*b(1:3))
            end do
        end if
    end do
end function create_grad_contra_ah2_operator

function create_grad_contra_ah_sbp_operator(domain, grad_operator_name) result(grad)

    use grad_contra_ah_sbp_mod,    only : grad_contra_ah_sbp_t
    use exchange_factory_mod,      only : create_symm_halo_exchange_Ah
    use ecs_metric_mod,            only : ecs_b1_proto, ecs_b2_proto, ecs_a1_proto, ecs_a2_proto
    use const_mod,                 only : pi

    type(domain_t),   intent(in)  :: domain
    character(len=*), intent(in)  :: grad_operator_name
    type(grad_contra_ah_sbp_t)    :: grad

    integer(kind=4)               :: halo_width_interior
    integer(kind=4), parameter    :: halo_width_edges=1

    integer(kind=4) :: t, is, ie, js, je, i, j
    real(kind=8)    :: alpha, beta, a(3), b(3)

    select case(grad_operator_name)
    case ("gradient_ah42_sbp_ecs")
        halo_width_interior = 3
        grad%sbp_operator_name="d42"
    case ("gradient_ah43_sbp_ecs")
        halo_width_interior = 5
        grad%sbp_operator_name="d43"
    case default
        call parcomm_global%abort("grad_factory_mod, create_grad_ah_sbp_operator"// &
                                  " - unknown SBP operator: "//grad_operator_name)
    end select

    grad%exch_scalar_interior =  &
              create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width_interior, 'full')
    grad%exch_vector_edges =  &
              create_symm_halo_exchange_Ah(domain%partition, domain%parcomm, &
                                         domain%topology,  halo_width_edges, 'full')

    grad%subtype = grad_operator_name

    allocate(grad%q(domain%mesh_xy%ts : domain%mesh_xy%te))
    do t=domain%mesh_xy%ts, domain%mesh_xy%te
        is = domain%mesh_xy%tile(t)%is
        ie = domain%mesh_xy%tile(t)%ie
        js = domain%mesh_xy%tile(t)%js
        je = domain%mesh_xy%tile(t)%je

        if(js == 1) then
            allocate(grad%q(t)%qb(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b1_proto(alpha,-0.25_8*pi)
                a(1:3) = ecs_a2_proto(alpha, 0.25_8*pi)
                a(1:3) = [a(1),-a(3),a(2)]
                grad%q(t)%qb(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(je == domain%mesh_xy%tile(t)%ny+1) then
            allocate(grad%q(t)%qt(is:ie))
            do i = is,ie
                alpha = -0.25_8*pi+(i-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b1_proto(alpha, 0.25_8*pi)
                a(1:3) = ecs_a2_proto(alpha,-0.25_8*pi)
                a(1:3) = [a(1),a(3),-a(2)]
                grad%q(t)%qt(i) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(is == 1) then
            allocate(grad%q(t)%ql(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b2_proto(-0.25_8*pi, beta)
                a(1:3) = ecs_a1_proto( 0.25_8*pi, beta)
                a(1:3) = [-a(3),a(2),a(1)]
                grad%q(t)%ql(j) = sum(a(1:3)*b(1:3))
            end do
        end if
        if(ie == domain%mesh_xy%tile(t)%nx+1) then
            allocate(grad%q(t)%qr(js:je))
            do j = js,je
                beta = -0.25_8*pi+(j-1)*domain%mesh_xy%tile(t)%hx
                b(1:3) = ecs_b2_proto( 0.25_8*pi, beta)
                a(1:3) = ecs_a1_proto(-0.25_8*pi, beta)
                a(1:3) = [a(3),a(2),-a(1)]
                grad%q(t)%qr(j) = sum(a(1:3)*b(1:3))
            end do
        end if
    end do

end function create_grad_contra_ah_sbp_operator

end module grad_factory_mod
