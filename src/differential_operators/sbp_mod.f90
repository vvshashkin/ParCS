module sbp_mod

use parcomm_mod,            only : parcomm_global

implicit none

private
public :: sbp_diff, sbp_apply

!Ah 2-1 scheme written as SBP
real(kind=8), parameter :: Q21(2,1) = reshape( [-1._8, 1._8],  [2,1])
integer(kind=4), parameter :: lastnonzeroQ21(1) =[2]
!2-th order diff non-staggered inner stencil and shift
real(kind=8),    parameter :: Da2_in(3) = [-0.5_8,0._8,0.5_8]
integer(kind=4), parameter :: Da2_inshift = -1

!boundary block of SBP diff matrix (inv(H)*Q)
real(kind=8), parameter :: Q42(6,4) = reshape( &
[-1.4117647058823528_8,   1.7352941176470589_8, -0.2352941176470588_8, -0.08823529411764705_8, 0.0_8,                 0.0_8, &
 -0.5_8,                  0.0_8,                0.5_8,                  0.0_8,                 0.0_8,                 0.0_8, &
  0.09302325581395347_8, -0.686046511627907_8,  0.0_8,                  0.686046511627907_8,  -0.09302325581395347_8, 0.0_8, &
  0.030612244897959186_8,-0.0_8,               -0.6020408163265307_8,   0.0_8,                 0.653061224489796_8,  -0.0816326530612245_8], &
  [6,4])
integer(kind=4), parameter :: lastnonzeroQ42(4) =[4,3,5,6]

real(kind=8), parameter :: Q43(6,4) = reshape( &
[-1.8280236842718398_8,   2.978054584697369_8,  -1.4653148294044123_8,  0.3078538227474199_8,  0.008136925288119662_8, -0.0007068190566565205_8, &
- 0.37814159126691327_8, -0.31480210500706657_8, 0.707290999364066_8,   0.04835554461933453_8,-0.06866771096803422_8,   0.005964863258613564_8,  &
  0.1162478575921292_8,  -0.8028500807354999_8,  0.21558841368737428_8, 0.5078566674295846_8, -0.03231754093993847_8,  -0.004525317033649745_8,  &
 -0.010287398573128293_8, 0.12592344801071614_8,-0.7341531396449149_8,  0.04979271660173103_8, 0.6506171865540596_8,  -0.08189281294846368_8],   &
                                              [6,4])
integer(kind=4), parameter :: lastnonzeroQ43(4) =[6,6,6,6]

!4-th order diff non-staggered inner stencil and shift
real(kind=8),    parameter :: Da4_in(5) = [1._8/12._8,-2._8/3.,0._8,2._8/3._8,-1._8/12._8]
integer(kind=4), parameter :: Da4_inshift = -2

contains

subroutine sbp_diff(df,is1,ie1,js1,je1,is,ie,js,je, &
                    direction,opername,f,is0,ie0,js0,je0,nsrc,ntarg)
    !dimensions of input function
    character(len=*), intent(in)    :: direction, opername
    integer(kind=4),  intent(in)    :: is0, ie0, js0, je0
    !global source and target mesh dimensions in the direction of differentiation:
    integer(kind=4),  intent(in)    :: nsrc, ntarg
    !output function dimensions:

    real(kind=8),     intent(in)    :: f(is0:ie0,js0:je0)

    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1
    integer(kind=4),  intent(in)    :: is,ie,js,je
    real(kind=8),     intent(inout) :: df(is1:ie1,js1:je1)

    if(direction /= 'x' .and. direction /= 'y') then
        call parcomm_global%abort("sbp_mod, sbp_diff - unknown differentiation direction: "// &
                                  direction // " . Use x or y")
    end if

    select case(opername)
    case ("d21")
        call sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je, &
                       f,is0,ie0,js0,je0,nsrc,ntarg,   &
                       Q21,lastnonzeroQ21,Da2_in,Da2_inshift,-1.0_8, &
                       direction)

    case ("d42")
        call sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je, &
                       f,is0,ie0,js0,je0,nsrc,ntarg,   &
                       Q42,lastnonzeroQ42,Da4_in,Da4_inshift,-1.0_8, &
                       direction)

    case ("d43")
        call sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je, &
                       f,is0,ie0,js0,je0,nsrc,ntarg,   &
                       Q43,lastnonzeroQ43,Da4_in,Da4_inshift,-1.0_8, &
                       direction)
    case default
        call parcomm_global%abort("sbp_mod, sbp_diff - unkonwn operator: "//opername)
    end select

end subroutine sbp_diff

subroutine sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je, &
                     f,is0,ie0,js0,je0,nsrc,ntarg,       &
                     Q,last_nonzero,Qin,inshift,right_edge_sign,direction)
    !dimensions of input function
    real(kind=8),     intent(in)    :: Q(:,:)
    integer(kind=4),  intent(in)    :: last_nonzero(:)
    real(kind=8),     intent(in)    :: Qin(:)
    integer(kind=4),  intent(in)    :: inshift
    integer(kind=4),  intent(in)    :: is0, ie0, js0, je0
    !global mesh dimension in x:
    integer(kind=4),  intent(in)    :: nsrc, ntarg
    !output function dimensions:
    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1
    integer(kind=4),  intent(in)    :: is, ie, js, je

    real(kind=8),     intent(in)    :: f(is0:ie0,js0:je0)
    character(len=*), intent(in)    :: direction
    real(kind=8),     intent(in)  :: right_edge_sign

    real(kind=8),     intent(inout) :: df(is1:ie1,js1:je1)

    integer(kind=4) :: i, j, mx, my, bs, inwidth, inend

    bs = size(Q,2)
    inwidth = size(Qin)
    inend = inwidth+inshift-1

    select case(direction)
    case('x')

    do j=js,je
        do i = max(is,1),min(bs,ie)
            mx = last_nonzero(i)
            df(i,j) = sum(Q(1:mx,i)*f(1:mx,j))
        end do
        do i = max(bs+1,is), min(ie,ntarg-bs)
            df(i,j) = sum(Qin(1:inwidth)*f(i+inshift:i+inend,j))
        end do
        do i=max(is,ntarg-bs+1),min(ie,ntarg)
            mx = last_nonzero(ntarg-i+1)
            df(i,j) = right_edge_sign*sum(Q(mx:1:-1,ntarg-i+1)*f(nsrc-mx+1:nsrc,j))
        end do
    end do

    case('y')

    do j = max(js,1),min(bs,je)
        do i=is,ie
            my = last_nonzero(j)
            df(i,j) = sum(Q(1:my,j)*f(i,1:my))
        end do
    end do
    do j = max(bs+1,js), min(je,ntarg-bs)
        do i=is,ie
            df(i,j) = sum(Qin(1:inwidth)*f(i,j+inshift:j+inend))
        end do
    end do
    do j=max(js,ntarg-bs+1),min(je,ntarg)
        do i=is,ie
            my = last_nonzero(ntarg-j+1)
            df(i,j) = right_edge_sign*sum(Q(my:1:-1,ntarg-j+1)*f(i,nsrc-my+1:nsrc))
        end do
    end do

    case default
    call parcomm_global%abort("unknown direction in sbp_apply: "//direction)
    end select
end subroutine sbp_apply

end module sbp_mod
