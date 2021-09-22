!Module for equiangular cubed sphere metric parameters
module ecs_metric_mod

use metric_mod,                only : metric_t
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use parcomm_mod,               only : parcomm_global

implicit none

type, extends(metric_t) :: ecs_metric_t
    class(cubed_sphere_topology_t), allocatable :: topology
contains
    procedure :: point_r => ecs_point_r
    procedure :: a1      => ecs_a1
    procedure :: a2      => ecs_a2
    procedure :: b1      => ecs_b1
    procedure :: b2      => ecs_b2
    procedure :: Q       => ecs_Q
    procedure :: QI      => ecs_QI
    procedure :: G       => ecs_G
end type ecs_metric_t

contains

function ecs_proto2realface(topology,rotation_matrix,panel_ind,r) result(r1)
    !transform prototype-face 3d Cartesian coords to real (spherical / elliptical) face
    class(cubed_sphere_topology_t), intent(in) :: topology
    real(kind=8),                   intent(in) :: rotation_matrix(3,3)
    integer(kind=4),                intent(in) :: panel_ind
    real(kind=8),                   intent(in) :: r(3)
    !output
    real(kind=8) :: r1(3)
    !local
    real(kind=8) :: rtemp(3)


    rtemp(1) = topology%ex(1,panel_ind)*r(1)+topology%ey(1,panel_ind)*r(2)- &
                                                    topology%n(1,panel_ind)*r(3)
    rtemp(2) = topology%ex(2,panel_ind)*r(1)+topology%ey(2,panel_ind)*r(2)- &
                                                    topology%n(2,panel_ind)*r(3)
    rtemp(3) = topology%ex(3,panel_ind)*r(1)+topology%ey(3,panel_ind)*r(2)- &
                                                    topology%n(3,panel_ind)*r(3)
    r1(1) = sum(rotation_matrix(1:3,1)*rtemp(1:3))
    r1(2) = sum(rotation_matrix(1:3,2)*rtemp(1:3))
    r1(3) = sum(rotation_matrix(1:3,3)*rtemp(1:3))
end function ecs_proto2realface

function ecs_point_r(this,panel_ind,alpha,beta) result(r)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: r(3)

    r = ecs_vector(this,panel_ind,"r",alpha,beta)
end function ecs_point_r

function ecs_a1(this,panel_ind,alpha,beta) result(a1)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a1(3)

    a1 = ecs_vector(this,panel_ind,"a1",alpha, beta)
end function ecs_a1

function ecs_a2(this,panel_ind,alpha,beta) result(a2)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: a2(3)

    a2 = ecs_vector(this,panel_ind,"a2",alpha,beta)
end function ecs_a2

function ecs_b1(this,panel_ind,alpha,beta) result(b1)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: b1(3)

    b1 = ecs_vector(this,panel_ind,"b1",alpha,beta)
end function ecs_b1

function ecs_b2(this,panel_ind,alpha,beta) result(b2)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha,beta
    real(kind=8)                    :: b2(3)

    b2 = ecs_vector(this,panel_ind,"b2",alpha,beta)
end function ecs_b2

function ecs_vector(this,panel_ind,vector_type,alpha,beta) result(r)
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    character(len=*),    intent(in) :: vector_type
    real(kind=8),        intent(in) :: alpha,beta
    real(kind=8)                    :: r(3)

    select case(vector_type)
    case("r")
        r = ecs_point_r_proto(alpha,beta)
    case("a1")
        r = ecs_a1_proto(alpha,beta)
    case("a2")
        r = ecs_a2_proto(alpha,beta)
    case("b1")
        r = ecs_b1_proto(alpha,beta)
    case("b2")
        r = ecs_b2_proto(alpha,beta)
    case default
        call parcomm_global%abort("unknown vector type " // vector_type //&
                                 " in ecs_metric_mod")
    end select

    r = ecs_proto2realface(this%topology,this%rotation_matrix,panel_ind,r)

end function ecs_vector

function ecs_point_r_proto(alpha,beta) result(r)
    real(kind=8), intent(in)    :: alpha, beta
    real(kind=8)                :: r(3)
    !locals:
    real(kind=8) d

    !grid point at proto cube face: (assumes x = tan(alpha), y = tan(beta), z = 1)
    d = sqrt(1._8+tan(alpha)**2+tan(beta)**2)
    !coordinates at prototype spherical face: vec(r)/||r||, r = (x, y, z)
    r(1) = tan(alpha)/d; r(2) = tan(beta)/d; r(3) = 1._8/d

end function ecs_point_r_proto

function ecs_a1_proto(alpha,beta) result(a1)
    !dr / d alpha vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: a1(3)
    !local
    real(kind=8) ta, tb

    ta = tan(alpha);   tb = tan(beta)
    !d (xyz)^T/ d alpha
    a1(1) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a1(2) = -ta*tb*(1d0+ta**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a1(3) = -ta*(1d0+ta**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

end function ecs_a1_proto

function ecs_a2_proto(alpha,beta) result(a2)
    !dr / d beta vector at prototype face for given alpha&beta
    real(kind=8)              :: a2(3)
    real(kind=8), intent(in)  :: alpha, beta
    !local
    real(kind=8) ta, tb

    ta = tan(alpha);   tb = tan(beta)
    !d (xyz)^T/ d beta
    a2(1) = -ta*tb*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a2(2) = (1d0+ta**2d0)*(1d0+tb**2d0)/(1d0+ta**2d0+tb**2d0)**(3d0/2d0)
    a2(3) = -tb*(1d0+tb**2d0) / (1d0+ta**2d0+tb**2d0)**(3d0/2d0)

end function ecs_a2_proto

function ecs_b1_proto(alpha,beta) result(b1)
    !contravariant alpha vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: b1(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)
    !d (xyz)^T/ d alpha
    b1(1) = (1._8+tb**2-(ta*tb)**2/(1._8+ta**2)) / sigm
    b1(2) = 0._8
    b1(3) =-ta*(1._8+tb**2/(1+ta**2))/sigm

end function ecs_b1_proto

function ecs_b2_proto(alpha,beta) result(b2)
    !contravariant beta vector at prototype face for given alpha&beta
    real(kind=8), intent(in)  :: alpha, beta
    real(kind=8)              :: b2(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)
    !d (xyz)^T/ d beta
    b2(1) = 0._8
    b2(2) = ((1+ta**2)-(ta*tb)**2/(1+tb**2))/sigm
    b2(3) =-tb*(ta**2/(1+tb**2)+1._8)/sigm

end function ecs_b2_proto

function ecs_Q(this,panel_ind,alpha,beta) result(Q)
    !Calculate metric tensor
    !Q = |a1*a1  a1*a2| = |Q(1) Q(2)|
    !    |a1*a2  a2*a2|   |Q(2) Q(3)|
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: Q(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    Q(1) = (1._8+ta**2)**2*(1._8+tb**2)/sigm**4
    Q(2) =-(1._8+ta**2)*(1._8+tb**2)*ta*tb/sigm**4
    Q(3) = (1._8+ta**2)*(1._8+tb**2)**2/sigm**4

end function ecs_Q

function ecs_QI(this,panel_ind, alpha, beta) result(QI)
    !Calculate inverse metric tensor
    !QI = inv |a1*a1  a1*a2| = |QI(1) QI(2)|
    !         |a1*a2  a2*a2|   |QI(2) QI(3)|
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha,beta
    real(kind=8)                    :: QI(3)
    !local
    real(kind=8) ta, tb, sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    Qi(1) = sigm**2 / (1._8+ta**2)
    Qi(2) = sigm**2*ta*tb / ((1+ta**2)*(1+tb**2))
    Qi(3) = sigm**2 / (1._8+tb**2)

end function ecs_QI

function ecs_G(this,panel_ind,alpha,beta) result(G)
    !Compute sqrt of metric tensor
    class(ecs_metric_t), intent(in) :: this
    integer(kind=4),     intent(in) :: panel_ind
    real(kind=8),        intent(in) :: alpha, beta
    real(kind=8)                    :: G
    !local
    real(kind=8) ta,tb,sigm

    ta = tan(alpha);   tb = tan(beta)
    sigm = sqrt(1._8+ta**2+tb**2)

    G = (1._8+ta**2)*(1._8+tb**2) / sigm**3
end function ecs_G

end module ecs_metric_mod
