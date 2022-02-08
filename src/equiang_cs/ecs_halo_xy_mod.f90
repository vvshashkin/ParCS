!Module to interpolate values to virtual points beyond face edge
!Edge interpolations are treated as one following case:
!    __________
!   / source  /
!  / /^b     /
! / +-->a   /
!/_________/     ^z
!|         |     |  /^y
!|  ^b     |     | /
!|  | targ |     |/
!|  +-->a  |     +------->x
!|_________|

!Edge and corners numeration:
!4__________3
!| edge2    |
!|e        e|
!|d        d|
!|g        g|
!|3        4|
!1----------2
!   edge1

module ecs_halo_xy_mod

use halo_mod,          only : halo_t
use exchange_halo_mod, only : exchange_t
use parcomm_mod,       only : parcomm_global

implicit none

type, extends(halo_t) :: ecs_halo_xy_t

    integer(kind=4)                       :: ts, te
    class(exchange_t),        allocatable :: exch_halo
    type(ecs_tile_halo_xy_t), allocatable :: tile(:)
    logical :: is_z_interfaces

    contains

    procedure :: get_halo_scalar => get_ecs_halo

end type

type ecs_tile_halo_xy_t
  integer(kind=4) n          !number of real grid points along cubed-sphere face edge
  integer(kind=4) halo_width !number of rows in halo-zone
  integer(kind=4) wsx, wex     !starting and ending index of interpolation weights x-edges
  integer(kind=4) wsy, wey     !starting and ending index of interpolation weights y-edges
  logical lhalo(4), lcorn(4) !flags to perform haloprocedures at specific edge and corner
  real(kind=8)   , allocatable :: wx(:,:,:) !interpolation weights for edges along x-axis
  integer(kind=4), allocatable :: indx(:,:) !interpolation stencil base point index for edges along x-axis
                                           !   s-----s-t----s------s
                                           !   base--^
  real(kind=8)   , allocatable :: wy(:,:,:) !interpolation weights, edges along y-axis
  integer(kind=4), allocatable :: indy(:,:) !interpolation stencil base point index, edgws along x-axis

  contains
  procedure, public :: interp  => ext_ecs_tile_halo

end type ecs_tile_halo_xy_t

contains


subroutine get_ecs_halo(this,f,domain,halo_width)
    use grid_field_mod, only : grid_field_t
    use domain_mod,     only : domain_t

    class(ecs_halo_xy_t),     intent(inout) :: this
    class(grid_field_t),      intent(inout) :: f
    type(domain_t),           intent(in)    :: domain
    integer(kind=4),          intent(in)    :: halo_width

    integer(kind=4) t

    call this%exch_halo%do(f, domain%parcomm)

    do t=this%ts,this%te
        if(this%is_z_interfaces) then
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_xyz%tile(t),halo_width)
        else
            call this%tile(t)%interp(f%tile(t),domain%partition%tiles_xy%tile(t),halo_width)
        end if
    end do
end subroutine get_ecs_halo

subroutine ext_ecs_tile_halo(this, f, tile, halo_width)
!interpolate source face values at halo zones to target face virtual points
use grid_field_mod, only : tile_field_t
use tile_mod,       only : tile_t

class(ecs_tile_halo_xy_t), intent(in)    :: this
type(tile_field_t),        intent(inout) :: f
type(tile_t),              intent(in)    :: tile
integer(kind=4),           intent(in)    :: halo_width

!local
integer(kind=4) n, thw, ihw
integer(kind=4) is, ie, js, je
integer(kind=4) nvi, nvj, nvk
integer(kind=4) isv, iev, jsv,jev, ksv, kev
integer(kind=4) klev
integer(kind=4) i,j,k
logical lhalo(4) !halo-procedure at edge

!short names for needed params
n = this%n
thw = this%halo_width !dimensions of initialized weights, theoretically max hw
ihw = halo_width      !width of user-requested halo-interpolations

if(thw < ihw) call parcomm_global%abort("requested halo width is greater than " // &
                                          "initialized maximum halo_width")

is = tile%is;      ie = tile%ie
js = tile%js;      je = tile%je
nvi = tile%is-f%is;    nvj = tile%js-f%js; nvk = tile%ks-f%ks

isv = is-nvi;   iev = ie+nvi
jsv = js-nvj;   jev = je+nvj
ksv = tile%ks-nvk; kev = tile%ke+nvk

klev = kev-ksv+1

lhalo = this%lhalo

if(lhalo(1).or.lhalo(2)) then
    call ecs_halo_edges_x(f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%wx,this%wsx,this%wex,thw,this%indx,      & !weight array and its bounds
                          lhalo,n,ihw)                                  !operating parameters
end if
if(lhalo(3).or.lhalo(4)) then
    call ecs_halo_edges_y(f%p,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                          this%wy,this%wsy,this%wey,thw,this%indy,      & !weight array and its bounds
                          lhalo,n,ihw)                                  !operating parameters
end if

end subroutine ext_ecs_tile_halo

subroutine ecs_halo_edges_x(pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            wx,wsx,wex,whw,indx,                         & !weights and indices  arrays and their bounds
                            lhalo,n,hw)                                    !global operating parameters

!in-output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev)
!input
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
!interpolation weights
integer(kind=4), intent(in) :: wsx,wex,whw
real(kind=8),    intent(in) :: wx(-1:2,wsx:wex,whw) !interpolation weights
integer(kind=4), intent(in) :: indx(wsx:wex,whw)
!operting parameters
logical,         intent(in) :: lhalo(4)
integer(kind=4), intent(in) :: n, hw
!local
integer(kind=4) k, j, i
integer(kind=4) ish, ieh, jsh, jeh
logical lfail_hw
logical lfail_halo_long
real(kind=8) zbufx(isv:iev, hw, klev)

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = nvj< hw
if(lfail_hw) call parcomm_global%abort("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data along edges to perform interpolations
ish = max(1,is-hw); ieh = min(n,ie+hw)
lfail_halo_long = .false.
lfail_halo_long = (minval(indx(ish,1:hw))-1<isv .or. maxval(indx(ieh,1:hw))+2>iev)
if(lfail_halo_long) call parcomm_global%abort("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points
!store them in sepparate arrays to not to spoil original f%p

!Halo-procedures along edges
if(lhalo(1)) then
    do k = 1,klev
        if(is==1) pf(1,1,k) = (pf(1,1,k)+pf(1,0,k)+pf(0,1,k)) / 3.0_8
        do i = max(is,2),min(ie,n-1)
            pf(i,1,k) = 0.5_8*(pf(i,1,k)+pf(i,0,k))
        end do
        if(ie==n) pf(n,1,k) = (pf(n,1,k)+pf(n,0,k)+pf(n+1,1,k)) / 3.0_8
    end do
    zbufx(isv:iev,1:hw,1:klev) = pf(isv:iev, -1:0-hw:-1,1:klev)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,wx,wsx,wex,whw,indx,n)
    pf(is:ie,1-hw:0,1:klev) = zbufx(is:ie,hw:1:-1,1:klev)
end if
if(lhalo(2)) then
    do k = 1,klev
        if(is == 1) pf(1,n,k) = (pf(1,n,k)+pf(1,n+1,k)+pf(0,n,k)) / 3.0_8
        do i = max(2,is),min(ie,n-1)
            pf(i,n,k) = 0.5_8*(pf(i,n,k)+pf(i,n+1,k))
        end do
        if(ie == n) pf(n,n,k) = (pf(n,n,k)+pf(n,n+1,k)+pf(n+1,n,k)) / 3.0_8
    end do
    zbufx(isv:iev,1:hw,1:klev) = pf(isv:iev, n+2:n+hw+1, 1:klev)
    call ecs_ext_halo_1e(zbufx,isv,iev,is,ie,hw,klev,wx,wsx,wex,whw,indx,n)
    pf(isv:iev, n+1:n+hw, 1:klev) = zbufx(isv:iev,1:hw,1:klev)
end if

end subroutine ecs_halo_edges_x

subroutine ecs_halo_edges_y(pf,is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj, & !function values array and its bounds
                            wy,wsy,wey,whw,indy,                         & !weights and indices  arrays and their bounds
                            lhalo,n,hw)                                    !global operating parameters

!in-output
real(kind=8), intent(inout) :: pf(isv:iev,jsv:jev,klev)
!input
integer(kind=4), intent(in) :: is,ie,js,je,isv,iev,jsv,jev,klev,nvi,nvj
!interpolation weights
integer(kind=4), intent(in) :: wsy,wey,whw
real(kind=8),    intent(in) :: wy(-1:2,wsy:wey,whw) !interpolation weights
integer(kind=4), intent(in) :: indy(wsy:wey,whw)
!operting parameters
logical,         intent(in) :: lhalo(4)
integer(kind=4), intent(in) :: n, hw
!local
integer(kind=4) k, j, i
integer(kind=4) jsh, jeh
logical lfail_hw
logical lfail_halo_long
real(kind=8) zbufy(jsv:jev, hw, klev)

!Internal checks to avoid (at least some part of) out-of-array errors
!check if we have enough data for edge halo-procedure
lfail_hw = nvj< hw
if(lfail_hw) call parcomm_global%abort("cubed sphere halo_width > grid_function_t halo width, can't continue")
!check if we have enough data along edges to perform interpolations
jsh = max(1,js-hw); jeh = min(n,je+hw)
lfail_halo_long = (minval(indy(jsh,1:hw))-1<jsv .or. maxval(indy(jeh,1:hw))+2>jev)
if(lfail_halo_long) call parcomm_global%abort("not enough points along cub.sph edge to perform halo-interpolations (nvi=5 needed)")

!calculate values at corner special points
!store them in sepparate arrays to not to spoil original f%p

!Halo-procedures along edges
if(lhalo(3)) then
    do k=1,klev
        do j=max(2,js),min(n-1,je)
            pf(1,j,k) = 0.5_8*(pf(1,j,k)+pf(0,j,k))
        end do
        do j=jsv, jev
            zbufy(j,1:hw,k) = pf(-1:0-hw:-1,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,wy,wsy,wey,whw,indy,n)
    do k=1,klev
        do j=1,hw
            pf(1-j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if
if(lhalo(4)) then
    do k=1, klev
        do j=max(2,js),min(je,n-1)
            pf(n,j,k) = 0.5_8*(pf(n,j,k)+pf(n+1,j,k))
        end do
        do j=jsv,jev
            zbufy(j,1:hw,k) = pf(n+2:n+hw+1,j,k)
        end do
    end do
    call ecs_ext_halo_1e(zbufy,jsv,jev,js,je,hw,klev,wy,wsy,wey,whw,indy,n)
    do k=1,klev
        do j=1,hw
            pf(n+j,jsv:jev,k) = zbufy(jsv:jev,j,k)
        end do
    end do
end if

end subroutine ecs_halo_edges_y

subroutine ecs_ext_halo_1e(zf,i1v,i2v,i1,i2,hw,klev,w,ws,we,whw,ind,n)
!input
integer(kind=4), intent(in) :: i1v,i2v,i1,i2,hw,klev
!in-output:
real(kind=8), intent(inout) :: zf(i1v:i2v,hw,klev)!input: source face values, output: interpolated target face values
!input
integer(kind=4), intent(in) :: ws,we,whw
real(kind=8),    intent(in) :: w(-1:2,ws:we,whw)
integer(kind=4), intent(in) :: ind(ws:we,whw)
integer(kind=4), intent(in) :: n
!locals
real(kind=8) zh(i1-hw:i2+hw)!buffer for interpolated values
integer i, j, k, ii

do k=1, klev
    do j=1, hw
        do i = i1,i2
            ii = ind(i,j)
            zh(i) = sum(w(:,i,j)*zf(ii-1:ii+2,j,k))
        end do
        zf(i1:i2,j,k) = zh(i1:i2)
    end do
end do

end subroutine ecs_ext_halo_1e

end module ecs_halo_xy_mod
