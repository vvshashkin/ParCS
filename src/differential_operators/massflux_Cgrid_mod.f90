module massflux_Cgrid_mod

use abstract_massflux_mod, only : massflux_operator_t
use grid_field_mod,        only : grid_field_t
use domain_mod,            only : domain_t
use halo_mod,              only : halo_t, halo_vec_t
use parcomm_mod,           only : parcomm_global

implicit none

type, extends(massflux_operator_t), public :: massflux_chalo_t

    integer(kind=4)                :: order
    class(halo_t),     allocatable :: halo
    class(halo_vec_t), allocatable :: halo_flux

    contains

    procedure :: calc_massflux => calc_c2_massflux

end type massflux_chalo_t

type, extends(massflux_operator_t), public :: massflux_c_sbp21_t

    contains

    procedure :: calc_massflux => calc_c_sbp21_massflux

end type massflux_c_sbp21_t

type, extends(massflux_operator_t), public :: massflux_c_sbp42_t

    contains

    procedure :: calc_massflux => calc_c_sbp42_massflux

end type massflux_c_sbp42_t

contains

subroutine calc_c2_massflux(this, fx, fy, f, u, v, domain)
    class(massflux_chalo_t), intent(inout) :: this
    type(domain_t),          intent(in)    :: domain
    type(grid_field_t),      intent(inout) :: f, u, v
    !output:
    type(grid_field_t),      intent(inout) :: fx, fy

    integer(kind=4) :: t


    select case(this%order)
    case(2)
        call this%halo%get_halo_scalar(f,domain,1)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c2_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case(4)
        call this%halo%get_halo_scalar(f,domain,1)
        do t = domain%mesh_p%ts, domain%mesh_p%te
            call calc_c4_massflux_tile(fx%tile(t), fy%tile(t), &
                                       f%tile(t), u%tile(t), v%tile(t), &
                                       domain%mesh_x%tile(t), domain%mesh_y%tile(t))
        end do
    case default
        call parcomm_global%abort("massflux_chalo_t is currently implemented only for orders=2,4")
    end select
    call this%halo_flux%get_halo_vector(fx,fy,domain,1)

end subroutine calc_c2_massflux

subroutine calc_c2_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i-1,j,k))
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i,j-1,k))
        end do; end do

    end do
end subroutine calc_c2_massflux_tile

subroutine calc_c4_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je; do i=is,ie
            fx%p(i,j,k) = u%p(i,j,k)*(-f%p(i-2,j,k)+7._8*f%p(i-1,j,k)+7._8*f%p(i,j,k)-f%p(i+1,j,k))/12._8
        end do; end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        do j=js,je; do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*(-f%p(i,j-2,k)+7._8*f%p(i,j-1,k)+7._8*f%p(i,j,k)-f%p(i,j+1,k))/12._8
        end do; end do

    end do
end subroutine calc_c4_massflux_tile

subroutine calc_c_sbp21_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_c_sbp21_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f, u, v
    !output:
    type(grid_field_t),        intent(inout) :: fx, fy

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%te

        call calc_c_sbp21_massflux_tile(fx%tile(t), fy%tile(t), &
                                        f%tile(t), u%tile(t), v%tile(t), &
                                        domain%mesh_x%tile(t), domain%mesh_y%tile(t),&
                                        domain%mesh_o%tile(t))

    end do

end subroutine calc_c_sbp21_massflux

subroutine calc_c_sbp21_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y, mesh_o)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke

    ks = mesh_x%ks; ke = mesh_x%ke

    do k=ks,ke
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        do j=js,je
            if(is == 1) fx%p(1,j,k) = u%p(1,j,k)*f%p(1,j,k)
            do i=max(is,2),min(ie,mesh_o%nx)
                fx%p(i,j,k) = u%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i-1,j,k))
            end do
            if(ie == mesh_o%nx+1) fx%p(mesh_o%nx+1,j,k) = u%p(mesh_o%nx+1,j,k)*f%p(mesh_o%nx,j,k)
        end do

        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        if(js == 1) then
            do i=is,ie
                fy%p(i,1,k) = v%p(i,1,k)*f%p(i,1,k)
            end do
        end if

        do j=max(js,2),min(je,mesh_o%ny); do i=is,ie
            fy%p(i,j,k) = v%p(i,j,k)*0.5_8*(f%p(i,j,k)+f%p(i,j-1,k))
        end do; end do

        if(je == mesh_o%ny+1) then
            do i=is,ie
                fy%p(i,mesh_o%ny+1,k) = v%p(i,mesh_o%ny+1,k)*f%p(i,mesh_o%ny,k)
            end do
        end if

    end do
end subroutine calc_c_sbp21_massflux_tile

subroutine calc_c_sbp42_massflux(this, fx, fy, f, u, v, domain)

    class(massflux_c_sbp42_t), intent(inout) :: this
    type(domain_t),            intent(in)    :: domain
    type(grid_field_t),        intent(inout) :: f, u, v
    !output:
    type(grid_field_t),        intent(inout) :: fx, fy

    integer(kind=4) :: t

    do t = domain%mesh_p%ts, domain%mesh_p%te

        call calc_c_sbp42_massflux_tile(fx%tile(t), fy%tile(t), &
                                        f%tile(t), u%tile(t), v%tile(t), &
                                        domain%mesh_x%tile(t), domain%mesh_y%tile(t),&
                                        domain%mesh_o%tile(t))

    end do

end subroutine calc_c_sbp42_massflux

subroutine calc_c_sbp42_massflux_tile(fx,fy,f,u,v,mesh_x,mesh_y, mesh_o)

    use sbp_mod, only : sbp_apply

    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t), intent(inout) :: fx, fy
    type(tile_field_t), intent(in)    :: f, u, v
    type(tile_mesh_t),  intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: k, ks, ke
    integer(kind=4) :: i, j
    integer(kind=4) :: is, ie, js, je
    integer(kind=4) :: isf, ief, jsf, jef
    integer(kind=4) :: isv, iev, jsv, jev

    real(kind=8)    :: fp(f%is:f%ie,f%js:f%je)

     !interpolation from p-points to vector points
    real(kind=8), parameter :: wp_edge(5,4) = reshape( &
    [0.9988019158947662_8,  0.37629049812372334_8, -0.24898674393171474_8, -0.12610567008674245_8,  0.0_8, &
     0.5893172370190968_8,  0.2849441265403871_8,   0.16216003586193575_8, -0.03642139942141965_8,  0.0_8, &
     0.21710064816314334_8, 0.23237013413978796_8,  0.3839577872309932_8,   0.16657143046607598_8,  0.0_8, &
    -0.18778023255417622_8, 0.17820763584084146_8,  0.6435451442907053_8,   0.42940773411277094_8, -0.06338028169014084_8],&
        [5,4])
    !optimized only vector to scalar interp
    ! real(kind=8), parameter :: wp_edge(5,4) = reshape( &
    ! [1.3630402181345285_8, -0.17804474904273135_8, -0.23303115631809243_8,  0.048035687226327554_8, 0.0_8, &
    !  0.46340770044238916_8, 0.4765661872892848_8,   0.15664452409426372_8, -0.09661841182593758_8,  0.0_8, &
    !  0.07545241951434666_8, 0.44794495248229815_8,  0.37775283649236274_8,  0.09884979151099299_8,  0.0_8, &
    ! -0.04413695843145285_8,-0.04040344754874627_8,  0.649837488701711_8,    0.4980831989686297_8,  -0.06338028169014084_8],&
    !    [5,4])
    !optimize both vector to scalar and scalar to vector interp
    ! real(kind=8), parameter :: wp_edge(5,4) = reshape( &
    ! [1.5033318188124725_8, -0.3057736326590643_8, -0.3984481911192596_8, 0.2008900049658833_8, 0.0_8, &
    !  0.43098284048835733_8, 0.5061952466249993_8, 0.194660985284931_8, -0.13183907239828724_8, 0.0_8, &
    ! -0.01526517971334979_8, 0.5302965998626628_8, 0.485202339414723_8, -0.00023375956403533802_8, 0.0_8, &
    !  0.029523830043030195_8, -0.1073451894687714_8, 0.5627386071183124_8, 0.5784630339975705_8, -0.06338028169014084_8],&
    !    [5,4])
    integer(kind=4), parameter :: wp_last_nonzero(4) = [4,4,4,5]
    !inner interpolation stencil
    real(kind=8), parameter :: w_in(4) = [-1._8/16._8, 9._8/16._8, 9._8/16._8, -1._8/16._8]

    do k = mesh_o%ks, mesh_o%ke

        is = mesh_o%is; ie = mesh_o%ie
        js = mesh_o%js; je = mesh_o%je

        do j=js,je
            do i=is,ie
                fp(i,j) = f%p(i,j,k)!mesh_o%G(i,j)*f%p(i,j,k)
            end do
        end do

        isv = fx%is; iev = fx%ie
        jsv = fx%js; jev = fx%je
        isf = f%is ; ief = f%ie
        jsf = f%js ; jef = f%je
        is = mesh_x%is; ie = mesh_x%ie
        js = mesh_x%js; je = mesh_x%je

        call sbp_apply(fx%p(isv:iev,jsv:jev,k),isv,iev,jsv,jev, &
                       is,ie,js,je,            &
                       fp(isf:ief,jsf:jef), isf, ief, jsf, jef, &
                       mesh_o%nx, mesh_o%nx+1,   &
                       wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'x')

        do j=js, je
            do i=is, ie
                fx%p(i,j,k) = u%p(i,j,k)*fx%p(i,j,k)! / mesh_x%G(i,j)
            end do
        end do

        isv = fx%is; iev = fx%ie
        jsv = fx%js; jev = fx%je
        isf = f%is ; ief = f%ie
        jsf = f%js ; jef = f%je
        is = mesh_y%is; ie = mesh_y%ie
        js = mesh_y%js; je = mesh_y%je

        call sbp_apply(fy%p(isv:iev,jsv:jev,k),isv,iev,jsv,jev, &
                       is,ie,js,je,            &
                       fp(isf:ief,jsf:jef), isf, ief, jsf, jef, &
                       mesh_o%ny, mesh_o%ny+1,   &
                       wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'y')

        do j=js, je
            do i=is, ie
                fy%p(i,j,k) = v%p(i,j,k)*fy%p(i,j,k)! / mesh_y%G(i,j)
            end do
        end do
    end do

end subroutine calc_c_sbp42_massflux_tile

end module massflux_Cgrid_mod
