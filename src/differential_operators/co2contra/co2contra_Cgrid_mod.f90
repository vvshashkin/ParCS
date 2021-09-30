module co2contra_Cgrid_mod

use abstract_co2contra_mod, only : co2contra_operator_t
use grid_field_mod,         only : grid_field_t, tile_field_t
use mesh_mod,               only : tile_mesh_t
use domain_mod,             only : domain_t
use exchange_abstract_mod,  only : exchange_t
use parcomm_mod,            only : parcomm_global

implicit none

type, extends(co2contra_operator_t), public :: co2contra_c_sbp_t
    character(len=:),  allocatable :: operator_name
    class(exchange_t), allocatable :: exchange_inner
    contains
        procedure :: transform => transform_co2contra_c_sbp
end type co2contra_c_sbp_t

contains

subroutine transform_co2contra_c_sbp(this, u_contra, v_contra, u_cov, v_cov, domain)
    class(co2contra_c_sbp_t), intent(inout) :: this
    type(domain_t),           intent(in)    :: domain
    type(grid_field_t),       intent(inout) :: u_cov, v_cov
    !output:
    type(grid_field_t),       intent(inout) :: u_contra, v_contra

    integer(kind=4) :: t

    call this%exchange_inner%do_vec(u_cov, v_cov, domain%parcomm)

    select case(this%operator_name)
    case ("co2contra_c_sbp21")
        do t=domain%mesh_o%ts, domain%mesh_o%te
            call transform_co2contra_c_sbp21_tile(u_contra%tile(t), v_contra%tile(t),&
                                                  u_cov%tile(t), v_cov%tile(t),      &
                                                  domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
                                                  domain%mesh_o%tile(t))
        end do
    case("co2contra_c_sbp42")
        do t=domain%mesh_o%ts, domain%mesh_o%te
            call transform_co2contra_c_sbp42_tile(u_contra%tile(t), v_contra%tile(t),&
                                                  u_cov%tile(t), v_cov%tile(t),      &
                                                  domain%mesh_u%tile(t), domain%mesh_v%tile(t), &
                                                  domain%mesh_o%tile(t))
        end do
    case default
        call parcomm_global%abort("unknown co2contra_c_sbp operator "// this%operator_name)
    end select
end subroutine transform_co2contra_c_sbp

subroutine transform_co2contra_c_sbp21_tile(u_contra, v_contra, u_cov, v_cov, mesh_u, mesh_v, mesh_o)

    type(tile_field_t), intent(in)    :: u_cov, v_cov
    type(tile_mesh_t),  intent(in)    :: mesh_u, mesh_v, mesh_o
    !output
    type(tile_field_t), intent(inout) :: u_contra, v_contra

    integer(kind=4) :: i, j, k, is, ie, js, je, ks, ke
    real(kind=8)    :: u_at_o(mesh_u%is:mesh_u%ie,mesh_u%js-1:mesh_u%je)
    real(kind=8)    :: v_at_o(mesh_u%is-1:mesh_u%ie,mesh_u%js:mesh_u%je)
    real(kind=8)    :: v_at_u, u_at_v

    ks = mesh_o%ks; ke = mesh_o%ke

    do k=ks, ke

        is = mesh_u%is; ie = mesh_u%ie
        js = mesh_u%js; je = mesh_u%je

        do j=js,je
            do i=max(is-1,1),min(ie,mesh_o%nx)
                v_at_o(i,j) = 0.5_8*(v_cov%p(i,j,k)+v_cov%p(i,j+1,k))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
            end do

            if(is == 1) then
                u_contra%p(1,j,k) = mesh_u%Qi(1,1,j)*u_cov%p(1,j,k)+v_at_o(1,j) / mesh_u%G(1,j)
            end if
            do i = max(is,2), min(ie,mesh_o%nx)
                v_at_u = 0.5_8*(v_at_o(i,j)+v_at_o(i-1,j))
                u_contra%p(i,j,k) = mesh_u%Qi(1,i,j)*u_cov%p(i,j,k)+v_at_u / mesh_u%G(i,j)
            end do
            if(ie == mesh_o%nx+1) then
                u_contra%p(ie,j,k) = mesh_u%Qi(1,ie,j)*u_cov%p(ie,j,k)+v_at_o(ie-1,j) / mesh_u%G(ie,j)
            end if
        end do

        is = mesh_v%is; ie = mesh_v%ie
        js = mesh_v%js; je = mesh_v%je

        do j=max(js-1,1),min(je,mesh_o%ny)
            do i = is,ie
                u_at_o(i,j) = 0.5_8*(u_cov%p(i+1,j,k)+u_cov%p(i,j,k))*mesh_o%G(i,j)*mesh_o%Qi(2,i,j)
            end do
        end do

        if(js == 1) then
            do i=is, ie
                v_contra%p(i,1,k) = mesh_v%Qi(3,i,1)*v_cov%p(i,1,k)+u_at_o(i,1) / mesh_v%G(i,1)
            end do
        end if
        do j = max(js,2), min(je,mesh_o%ny)
            do i = is, ie
                u_at_v = 0.5_8*(u_at_o(i,j)+u_at_o(i,j-1))
                v_contra%p(i,j,k) = mesh_v%Qi(3,i,j)*v_cov%p(i,j,k)+u_at_v / mesh_v%G(i,j)
            end do
        end do
        if(je == mesh_o%ny+1) then
            do i = is, ie
                v_contra%p(i,je,k) = mesh_v%Qi(3,i,je)*v_cov%p(i,je,k)+u_at_o(i,je-1) / mesh_v%G(i,je)
            end do
        end if

    end do

end subroutine transform_co2contra_c_sbp21_tile

subroutine transform_co2contra_c_sbp42_tile(gx, gy, dx, dy, mesh_x, mesh_y, mesh_o)

    use sbp_mod, only : sbp_apply

    type(tile_field_t),     intent(inout) :: gx, gy
    type(tile_field_t),     intent(inout)    :: dx, dy
    type(tile_mesh_t),      intent(in)    :: mesh_x, mesh_y, mesh_o

    integer(kind=4) :: i, j, k, l, nsrc, ntarg
    integer(kind=4) :: is, ie, js, je, ks, ke
    !indices for grad components:
    integer(kind=4) :: isx, iex, jsx, jex
    integer(kind=4) :: isy, iey, jsy, jey
    !indices for grad components interpolated to p-points:
    integer(kind=4) :: ispx, iepx, jspx, jepx
    integer(kind=4) :: ispy, iepy, jspy, jepy

    integer(kind=4), parameter :: hw = 3 !maximum of interpolation stencil displacement
    real(kind=8)    :: dx_at_p(mesh_y%is:mesh_y%ie,mesh_y%js-hw:mesh_y%je+hw)
    real(kind=8)    :: dy_at_p(mesh_x%is-hw:mesh_x%ie+hw,mesh_x%js:mesh_x%je)

    !Interpolation weights
    !interpolation from vector points to p-points
    !non-optimized option (minimum l2 of coefficients)
    !real(kind=8), parameter :: wvec_edge(6,4) = reshape( &
    ![0.3585442775006853_8,  0.6119832845967544_8,  0.20040059830444001_8, -0.1709281604018784_8,  0.0_8,    0.0_8,    &
    ! 0.16724022138832148_8, 0.36635673412335484_8, 0.2655658675883291_8,   0.20083717689999595_8, 0.0_8,    0.0_8,    &
    !-0.09295505106784016_8, 0.1751328387308906_8,  0.36859947574175345_8,  0.609222736595201_8,  -0.06_8,   0.0_8,    &
    !-0.04904109392262206_8,-0.04097407434909711_8, 0.16657143046607598_8,  0.42344373780564915_8, 0.5625_8,-0.0625_8],&
    !   [6,4])
    !Optimized vector->scalar x^2 interpolation
    real(kind=8), parameter :: wvec_edge(6,4) = reshape( &
    [ 0.4892964885611128_8,    0.48123107353632727_8, 0.0696483872440123_8, -0.040175949341450676_8, 0.0_8,     0.0_8, &
    -0.07913099957454726_8,   0.6127279550862234_8,  0.5119370885511979_8, -0.045534044062872786_8, 0.0_8,     0.0_8, &
    -0.0869982983587545_8,    0.1691760860218048_8,  0.3626427230326682_8,  0.6151794893042863_8,  -0.06_8,    0.0_8, &
    0.018680545032460714_8, -0.10869571330417979_8, 0.09884979151099299_8, 0.4911653767607321_8,   0.5625_8, -0.0625_8],&
        [6,4])
    !Optimized vector <-> scalar interp in x^2
    ! real(kind=8), parameter :: wvec_edge(6,4) = reshape( &
    !  [ 0.5396575759839645_8, 0.4475591035840635_8, -0.014090935120015191_8, 0.026874255551989024_8, 0.0_8, 0.0_8, &
    !   -0.13589939229291748_8, 0.6508224599464276_8, 0.6060532569859003_8, -0.12097632463940905_8, 0.0_8, 0.0_8,   &
    !   -0.14875399135119022_8, 0.2102338641077255_8, 0.46579424583813406_8, 0.5327258814053357_8, -0.06_8, 0.0_8,  &
    !    0.07812389082006573_8, -0.14831895644807316_8, -0.00023375956403533802_8, 0.5704288251920487_8, 0.5625_8, -0.0625_8],&
    !     [6,4])
    integer(kind=4), parameter :: wvec_last_nonzero(4) = [4,4,5,6]
    !interpolation from p-points to vector points
    ! real(kind=8), parameter :: wp_edge(5,4) = reshape( &
    ! [0.9988019158947662_8,  0.37629049812372334_8, -0.24898674393171474_8, -0.12610567008674245_8,  0.0_8, &
    !  0.5893172370190968_8,  0.2849441265403871_8,   0.16216003586193575_8, -0.03642139942141965_8,  0.0_8, &
    !  0.21710064816314334_8, 0.23237013413978796_8,  0.3839577872309932_8,   0.16657143046607598_8,  0.0_8, &
    ! -0.18778023255417622_8, 0.17820763584084146_8,  0.6435451442907053_8,   0.42940773411277094_8, -0.06338028169014084_8],&
    !    [5,4])
    !optimized only vector to scalar interp
    real(kind=8), parameter :: wp_edge(5,4) = reshape( &
    [1.3630402181345285_8, -0.17804474904273135_8, -0.23303115631809243_8,  0.048035687226327554_8, 0.0_8, &
    0.46340770044238916_8, 0.4765661872892848_8,   0.15664452409426372_8, -0.09661841182593758_8,  0.0_8, &
    0.07545241951434666_8, 0.44794495248229815_8,  0.37775283649236274_8,  0.09884979151099299_8,  0.0_8, &
    -0.04413695843145285_8,-0.04040344754874627_8,  0.649837488701711_8,    0.4980831989686297_8,  -0.06338028169014084_8],&
        [5,4])
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

    ks = mesh_o%ks; ke = mesh_o%ke

    !bounding box for required values of gx
    isx = mesh_x%is; iex = mesh_x%ie
    jsx = mesh_x%js; jex = mesh_x%je
    !bounding box for needed values of dy_at_p:
    ispy = interp_p_stencil_start(isx,mesh_o%nx)
    iepy = interp_p_stencil_end  (iex,mesh_o%nx)
    jspy = jsx
    jepy = jex

    !bounding box for required values of gy
    isy = mesh_y%is; iey = mesh_y%ie
    jsy = mesh_y%js; jey = mesh_y%je
    !bounding box for needed values of dx_at_p:
    ispx = isy
    iepx = iey
    jspx = interp_p_stencil_start(jsx,mesh_o%ny)
    jepx = interp_p_stencil_end  (jex,mesh_o%ny)

    do k = ks,ke

       !interpolate dy to p-points
       ntarg = mesh_o%ny
       nsrc  = mesh_o%ny+1
        call sbp_apply(dy_at_p,isx-hw,iex+hw,jsx,jex,ispy,iepy,jspy,jepy,           &
                      dy%p(dy%is:dy%ie,dy%js:dy%je,k), dy%is, dy%ie, dy%js, dy%je, &
                      mesh_o%ny+1,mesh_o%ny,   &
                      wvec_edge,wvec_last_nonzero,w_in,-1,1.0_8,'y')
        !mutiplicate by metric terms in p-points
        do j = jspy,jepy
            do i = ispy, iepy
                dy_at_p(i,j) = mesh_o%G(i,j)*mesh_o%Qi(2,i,j)*dy_at_p(i,j)
            end do
        end do
        !interpolate to u-points
        call sbp_apply(gx%p(gx%is:gx%ie,gx%js:gx%je,k),gx%is,gx%ie,gx%js,gx%je, &
                      isx,iex,jsx,jex,            &
                      dy_at_p, isx-hw,iex+hw,jsx,jex, &
                      mesh_o%nx,mesh_o%nx+1,   &
                      wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'x')
        !interpolate dx to p-points
        call sbp_apply(dx_at_p,isy,iey,jsx-hw,jex+hw,ispx,iepx,jspx,jepx,           &
                      dx%p(dx%is:dx%ie,dx%js:dx%je,k), dx%is, dx%ie, dx%js, dx%je, &
                      mesh_o%nx+1,mesh_o%nx,   &
                      wvec_edge,wvec_last_nonzero,w_in,-1,1.0_8,'x')
        !multiplicate by metric terms in p-points
        do j = jspx,jepx
            do i = ispx, iepx
                dx_at_p(i,j) = mesh_o%G(i,j)*mesh_o%Qi(2,i,j)*dx_at_p(i,j)
            end do
        end do
        !interpolate dx to v-points
        call sbp_apply(gy%p(gy%is:gy%ie,gy%js:gy%je,k),gy%is,gy%ie,gy%js,gy%je, &
                      isy,iey,jsy,jey,            &
                      dx_at_p, isy,iey,jsy-hw,jey+hw, &
                      mesh_o%ny,mesh_o%ny+1,   &
                      wp_edge,wp_last_nonzero,w_in,-2,1.0_8,'y')

        do j=jsx,jex
            do i=isx,iex
                gx%p(i,j,k) = mesh_x%Qi(1,i,j)*dx%p(i,j,k)+gx%p(i,j,k)/mesh_x%G(i,j)
            end do
        end do

        do j=jsy,jey
            do i=isy,iey
                gy%p(i,j,k) = mesh_y%Qi(3,i,j)*dy%p(i,j,k)+gy%p(i,j,k)/mesh_y%G(i,j)
            end do
        end do
    end do

    contains
    !from u/v points to p points
    integer(kind=4) function interp_v_stencil_start(i, n) result(is)
        integer(kind=4) :: i, n
        if(i <= 4) then
            is = 1
        else
            is = min(i-1,n-3)
        end if
    end

    integer(kind=4) function interp_v_stencil_end(i, n) result(ie)
        integer(kind=4) :: i, n
        if(i <= n-5) then
            ie = max(4,i+2)
        else
            ie = n
        end if
    end

    !from p-points to u/v points
    integer(kind=4) function interp_p_stencil_start(i, n) result(is)
        integer(kind=4) :: i, n
        if(i <= 4) then
            is = 1
        else
            is = min(i-2,n-3)
        end if
    end

    integer(kind=4) function interp_p_stencil_end(i, n) result(ie)
        integer(kind=4) :: i, n
        if(i <= n-3) then
            ie = max(4,i+1)
        else
            ie = n
        end if
    end
end subroutine transform_co2contra_c_sbp42_tile


end module co2contra_Cgrid_mod
