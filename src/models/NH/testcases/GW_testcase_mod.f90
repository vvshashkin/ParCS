module NH_GW_testcase_mod

use parcomm_mod,        only : parcomm_global
use domain_mod,         only : domain_t
use stvec_mod,          only : stvec_t
use stvec_nh_mod,       only : stvec_nh_t
use test_fields_3d_mod, only : scalar_field3d_t

implicit none

type, extends(scalar_field3d_t) :: GW_theta_t
    real(kind=8) :: amp = 0.01_8
    real(kind=8) :: a = 5e3, Lz = 10e3
    real(kind=8) :: hor_scale
    contains
    procedure :: get_scalar_field_tile => get_GWtheta_tile
end type GW_theta_t

contains

subroutine get_GW_initial_conditions(stvec, domain, testcase_name)
    class(stvec_t),   intent(inout) :: stvec
    type(domain_t),   intent(in)    :: domain
    character(len=*), intent(in)    :: testcase_name

    logical :: is_nonlin

    select case(testcase_name)
    case ("GW")
        is_nonlin = .true.
    case("GW_linear")
        is_nonlin = .false.
    case default
        call parcomm_global%abort("unknown type of NH GW test: "// testcase_name)
    end select

    select type(stvec)
    type is (stvec_nh_t)
        call get_GW_initial_conditions_nh_stvec(stvec, domain, is_nonlin)
    class default
        call parcomm_global%abort("NH GW testcases does not support this type of stvec")
    end select

end subroutine get_GW_initial_conditions

subroutine get_GW_initial_conditions_nh_stvec(stvec, domain, is_nonlin)

    use grid_field_mod,                only: grid_field_t
    use grid_field_factory_mod,        only: create_grid_field
    use vertical_test_field_mod,       only: vertical_ExnerP_t, vertical_theta_t
    use const_N_profile_mod,           only: const_N_profile_t

    class(stvec_nh_t),   intent(inout) :: stvec
    type(domain_t),      intent(in)    :: domain
    logical,             intent(in)    :: is_nonlin

    real(kind=8), parameter :: Nb = 0.01_8, T0 = 300.0_8, pref = 1e5_8

    type(GW_theta_t) :: theta_generator
    type(grid_field_t) :: P0, theta0
    type(vertical_ExnerP_t) :: P0_generator
    type(vertical_theta_t)  :: theta0_generator

    call stvec%u%assign(0.0_8,domain%mesh_u)
    call stvec%v%assign(0.0_8,domain%mesh_v)
    call stvec%eta_dot%assign(0.0_8,domain%mesh_w)
    call stvec%P%assign(0.0_8,domain%mesh_p)

    theta_generator = GW_theta_t(a=5e3,Lz=10e3,hor_scale=domain%mesh_w%scale)
    call theta_generator%get_scalar_field(stvec%theta,domain%mesh_w,0)

    if(is_nonlin) then
        call create_grid_field(P0,0,0,domain%mesh_p)
        call create_grid_field(theta0,0,0,domain%mesh_w)
        P0_generator%t0 = T0
        P0_generator%p0 = pref
        P0_generator%vert_profile = const_N_profile_t(N=Nb)
        theta0_generator%t0 = T0
        theta0_generator%p0 = pref
        theta0_generator%vert_profile = const_N_profile_t(N=Nb)

        call P0_generator%get_scalar_field(P0,domain%mesh_p,0)
        call theta0_generator%get_scalar_field(theta0,domain%mesh_w,0)

        call stvec%P%update(1.0_8,P0,domain%mesh_p)
        call stvec%theta%update(1.0_8,theta0,domain%mesh_w)
    end if

end subroutine get_GW_initial_conditions_nh_stvec

subroutine get_GWtheta_tile(this,f,mesh,halo_width)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t
    use const_mod,      only : pi

    class(GW_theta_t),    intent(in)    :: this
    type(tile_field_t),   intent(inout) :: f
    type(tile_mesh_t),    intent(in)    :: mesh
    integer(kind=4),      intent(in)    :: halo_width

    integer(kind=4) :: i,j,k,is,ie,js,je,ks,ke
    real(kind=8)    :: d

    is = mesh%is; ie = mesh%ie
    js = mesh%js; je = mesh%je
    ks = mesh%ks; ke = mesh%ke

    do k=ks,ke
        do j=js,je
            do i=is,ie
                d = acos(mesh%rx(i,j,k))*this%hor_scale
                f%p(i,j,k) = this%amp*this%a**2 / (this%a**2+d**2)*sin(pi*mesh%h(i,j,k)/this%Lz)
            end do
        end do
    end do

end subroutine get_GWtheta_tile

end module NH_GW_testcase_mod
