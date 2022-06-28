module test_metric_mod

implicit none

private
public :: test_metric, test_metric_vert

real(kind=8), parameter :: test_tolerance=3e-15

contains

subroutine test_metric()

use mpi

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use topology_mod,           only : topology_t
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use mesh_factory_mod,       only : create_mesh
use mesh_mod,               only : mesh_t
use config_domain_mod,      only : config_domain_t
use const_mod,              only : Earth_radii

type(domain_t)                     :: domain
type(grid_field_t)                 :: f1, f2
type(mesh_t),          allocatable :: mesh(:)
type(config_domain_t)              :: config_domain

integer(kind=4), parameter         :: nh=16, nz=10, halo_width=3
real(kind=8),    parameter         :: zeps = 1e-16_8
integer(kind=4)                    :: myid, np, ierr, code

integer(kind=4) :: ts, te
integer(kind=4) :: ind, i, j, k, ifc, err_sum, gl_err_sum
integer(kind=4) :: ie, is, js, je

integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number

character(:), allocatable :: frmt
logical :: lcross_edge_xyz = .true., lsymm_check = .true.
real(kind=8) zq(1:nh,3), zp(1:nh,3)
real(kind=8) err_ort_rab, err_ort_abt, err_q, err_qi, err_g

config_domain%N  = nh
config_domain%Nz = nz
config_domain%staggering_type     = "C"
config_domain%vertical_staggering = "CharneyPhilips"
config_domain%metric_type         = "shallow_atmosphere_metric"
config_domain%topology_type       = "cube"
config_domain%h_top = 30e3
call config_domain%config_metric%set_defaults()
config_domain%config_metric%scale = Earth_radii
config_domain%config_metric%vertical_scale = config_domain%h_top
config_domain%is_orographic_curvilinear = .true.
config_domain%orography_name = "test_orography"

call create_domain(domain, config_domain)
call domain%parcomm%print('equiangular cubed-sphere test')

if(domain%parcomm%np >1) then
    print *, "this test currently runs only with 1 mpi-process. exiting"
    return
end if

do ind = 1,1!domain%partition%ts, domain%partition%te
    lsymm_check = lsymm_check .and. &
        symmetricity_check_h(domain%topology, domain%mesh_p, domain%partition%tiles_p%tile(ind), &
                         ind, domain%partition%panel_map(ind)) .and. &
        symmetricity_check_h(domain%topology, domain%mesh_u, domain%partition%tiles_u%tile(ind), &
                         ind, domain%partition%panel_map(ind)) .and. &
        symmetricity_check_h(domain%topology, domain%mesh_v, domain%partition%tiles_v%tile(ind), &
                         ind, domain%partition%panel_map(ind))
end do

if(lsymm_check) then
    print *, "symmetricity test passed"
else
    print *, "symmetricity test failed"
end if

call test_mesh_metric(domain%mesh_p,halo_width,"p")
call test_mesh_metric(domain%mesh_u,halo_width,"u")
call test_mesh_metric(domain%mesh_v,halo_width,"v")

contains

end subroutine test_metric

!check symmetricity of points on ecs grid face
logical function symmetricity_check_h(topology,mesh,tile,ind,panel_ind) result(lsym_check)
    use topology_mod, only : topology_t
    use mesh_mod,     only : mesh_t
    use tile_mod,     only : tile_t

    class(topology_t), intent(in)      :: topology
    class(mesh_t),     intent(in)      :: mesh
    type(tile_t),      intent(in)      :: tile
    integer(kind=4), intent(in)        :: ind, panel_ind
    !local
    real(kind=8) ee(3)
    integer(kind=4) is, ie, js, je, ks, ke


    lsym_check = .true.

    call tile%getind(is,ie,js,je, ks, ke)

    ee = real(topology%ex(:,panel_ind),8)
    lsym_check = lsym_check .and. symmetric( is,ie,js,js,ks,ke,ee)
    lsym_check = lsym_check .and. symmetric( is,ie,je,je,ks,ke,ee)
    ee = real(topology%ey(:,panel_ind),8)
    lsym_check = lsym_check .and. symmetric(is,is,js,je,ks,ke,ee)
    lsym_check = lsym_check .and. symmetric(ie,ie,js,je,ks,ke,ee)

    contains
    logical function symmetric(i1,i2,j1,j2,k1,k2,ee) result(lsym)
        integer(kind=4), intent(in) :: i1,i2,j1,j2,k1,k2
        real(kind=8),    intent(in) :: ee(3)
        !local
        real(kind=8) a(3), b(3), zpr
        real(kind=8), parameter :: zeps = 1e-16_8

        a = [mesh%tile(ind)%rx(i1,j1,k1), mesh%tile(ind)%ry(i1,j1,k1), mesh%tile(ind)%rz(i1,j1,k1)]
        b = [mesh%tile(ind)%rx(i2,j2,k1), mesh%tile(ind)%ry(i2,j2,k1), mesh%tile(ind)%rz(i2,j2,k1)]
        zpr = sum(a*ee)
        a = a - 2._8*zpr*ee !flip a-vector
        lsym = sum(abs(a-b))<zeps
    end function symmetric

end function symmetricity_check_h

!checks the following facts
!1)cubsph basis vectors tangent to sphere
!2)covariant (a) vectors are orthogonal to contravariant (b) vectors
!3)metric tensor is matrix of covariant vectors dot-prods
!4)b1 = qi11*a1+qi12*a2 and the same for b2
!5)J = sqrt(det(Q)) is sqrt of determinant of Q
subroutine test_mesh_metric(mesh,halo_width,points_type)
    use mesh_mod, only : mesh_t

    class(mesh_t),    intent(in) :: mesh
    integer(kind=4),  intent(in) :: halo_width
    character(len=*), intent(in) :: points_type

    integer(kind=4) ts,te,is,ie,js,je,ks,ke
    integer(kind=4) i, j, k, ind
    real(kind=8) :: err_ort_rab=0.0_8, err_ort_abt=0.0_8, err_q=0.0_8, &
                    err_qi =0.0_8, err_g=0.0_8
    real(kind=8) :: r(3), a1(4), a2(4), a3(4), b1(3), b2(3), b3(4), q(6), qi(6), Jac

    ts = mesh%ts
    te = mesh%te

    do ind=ts,te
        is = mesh%tile(ind)%is
        ie = mesh%tile(ind)%ie
        js = mesh%tile(ind)%js
        je = mesh%tile(ind)%je
        ks = mesh%tile(ind)%ks
        ke = mesh%tile(ind)%ke

        do k = ks, ke
        do j=js-halo_width,je+halo_width
            do i=is-halo_width,ie+halo_width
                r(1) = mesh%tile(ind)%rx(i,j,k); r(2) = mesh%tile(ind)%ry(i,j,k); r(3) = mesh%tile(ind)%rz(i,j,k)
                a1(1:4) = mesh%tile(ind)%a1(1:4,i,j,k)
                a2(1:4) = mesh%tile(ind)%a2(1:4,i,j,k)
                a3(1:4) = mesh%tile(ind)%a3(1:4,i,j,k)
                b1(1:3) = mesh%tile(ind)%b1(1:3,i,j,k)
                b2(1:3) = mesh%tile(ind)%b2(1:3,i,j,k)
                b3(1:4) = mesh%tile(ind)%b3(1:4,i,j,k)
                q(1:6)  = mesh%tile(ind)%Q(1:6,i,j,k)
                qi(1:6) = mesh%tile(ind)%Qi(1:6,i,j,k)
                Jac     = mesh%tile(ind)%J(i,j,k)

                err_ort_rab = max(err_ort_rab, abs(sum(a1(1:3)*r(1:3))))
                err_ort_rab = max(err_ort_rab, abs(sum(a2(1:3)*r(1:3))))
                err_ort_rab = max(err_ort_rab, abs(sum(b1(1:3)*r(1:3))))
                err_ort_rab = max(err_ort_rab, abs(sum(b2(1:3)*r(1:3))))

                err_ort_abt = max(err_ort_abt, abs(1.0_8-sum(a1(1:3)*b1(1:3))))
                err_ort_abt = max(err_ort_abt, abs(sum(a1(1:3)*b2(1:3))))
                err_ort_abt = max(err_ort_abt, abs(sum(a1(1:4)*b3(1:4))))
                err_ort_abt = max(err_ort_abt, abs(1.0_8-sum(a2(1:3)*b2(1:3))))
                err_ort_abt = max(err_ort_abt, abs(sum(a2(1:3)*b1(1:3))))
                err_ort_abt = max(err_ort_abt, abs(sum(a2(1:4)*b3(1:4))))
                err_ort_abt = max(err_ort_abt, abs(1.0_8-sum(a3(1:4)*b3(1:4))))
                err_ort_abt = max(err_ort_abt, abs(sum(a3(1:3)*b1(1:3))))
                err_ort_abt = max(err_ort_abt, abs(sum(a3(1:3)*b2(1:3))))

                err_q = max(err_q, abs(sum(a1(1:4)**2)-q(1))+&
                                   abs(sum(a1(1:4)*a2(1:4))- q(2))+&
                                   abs(sum(a2(1:4)**2)-q(3)))

                err_qi = max(err_qi, sum(abs(b1(1:3)-qi(1)*a1(1:3)-qi(2)*a2(1:3))))
                err_qi = max(err_qi, sum(abs(b2(1:3)-qi(2)*a1(1:3)-qi(3)*a2(1:3))))
                err_qi = max(err_qi, abs(qi(1)*a1(4)+qi(2)*a2(4)+qi(4)*a3(4)))
                err_qi = max(err_qi, abs(qi(2)*a1(4)+qi(3)*a2(4)+qi(5)*a3(4)))

                err_g = max(err_g, abs(sqrt(q(1)*q(3)*q(6)+2.0_8*q(2)*q(4)*q(5)-     &
                                            q(1)*q(5)**2-q(2)**2*q(6)-q(3)*q(4)**2)- &
                                            Jac))
            end do
        end do
        end do
    end do
    print *, "metric errors @ "//points_type//"-points", err_ort_rab, err_ort_abt, err_q, err_qi, err_g
    if(max(err_ort_rab, err_ort_abt, err_q, err_qi, err_g) < test_tolerance) then
        print *, "shallow-atm ecs mesh-metric test passed ("//points_type//"-points)"
    else
        print *, "shallow-atm ecs mesh-metric test failed ("//points_type//"-points)"
    end if

end

subroutine test_metric_vert()
    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use config_domain_mod,      only : config_domain_t

    type(domain_t)        :: domain
    type(config_domain_t) :: config_domain
    integer(kind=4), parameter         :: nh=16, nz=3
    logical :: is_passed

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = "C"
    config_domain%vertical_staggering = "CharneyPhilips"
    config_domain%metric_type         = "shallow_atmosphere_metric"
    config_domain%topology_type       = "cube"
    config_domain%h_top = 30e3
    call config_domain%config_metric%set_defaults()

    call create_domain(domain, config_domain)

    is_passed = (domain%mesh_z%tile(1)%ke == domain%mesh_o%tile(1)%ke+1)
    is_passed = is_passed .and. (domain%mesh_z%tile(1)%hz == domain%mesh_o%tile(1)%hz)
    is_passed = is_passed .and. (domain%mesh_xyz%tile(1)%ke == domain%mesh_xy%tile(1)%ke+1)
    is_passed = is_passed .and. (domain%mesh_xyz%tile(1)%hz == domain%mesh_xy%tile(1)%hz)

    if(is_passed) then
        print *, "vertical metric test passed"
    else
        print *, "vertical metric test failed"
    end if
    ! print *, domain%mesh_z%tile(1)%h(1,1,:)
    ! print *, domain%mesh_xyz%tile(1)%h(1,1,:)

end subroutine test_metric_vert
end module test_metric_mod
