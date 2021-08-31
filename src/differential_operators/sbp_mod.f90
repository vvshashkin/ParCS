module sbp_mod

use parcomm_mod,            only : parcomm_global

implicit none

private
public :: sbp_diff

!boundary block of SBP diff matrix (inv(H)*Q)
real(kind=8), parameter :: Q42(6,4) = reshape( &
[-1.4117647058823528_8,   1.7352941176470589_8, -0.2352941176470588_8, -0.08823529411764705_8, 0.0_8,                 0.0_8, &
 -0.5_8,                  0.0_8,                0.5_8,                  0.0_8,                 0.0_8,                 0.0_8, &
  0.09302325581395347_8, -0.686046511627907_8,  0.0_8,                  0.686046511627907_8,  -0.09302325581395347_8, 0.0_8, &
  0.030612244897959186_8,-0.0_8,               -0.6020408163265307_8,   0.0_8,                 0.653061224489796_8,  -0.0816326530612245_8], &
  [6,4])
integer(kind=4), parameter :: lastnonzeroQ42(4) =[4,3,5,6]

contains

subroutine sbp_diff(direction,opername,f,is,ie,js,je,n,is1,ie1,js1,je1,df)
    !dimensions of input function
    character(len=*), intent(in)    :: direction, opername
    integer(kind=4),  intent(in)    :: is, ie, js, je
    !global mesh dimension in the direction of differentiation:
    integer(kind=4),  intent(in)    :: n
    !output function dimensions:
    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1

    real(kind=8),     intent(in)    :: f(is:ie,js:je)
    real(kind=8),     intent(inout) :: df(is1:ie1,js1:je1)

    if(direction /= 'x' .and. direction /= 'y') then
        call parcomm_global%abort("sbp_mod, sbp_diff - unknown differentiation direction: "// &
                                  direction // " . Use x or y")
    end if

    select case(opername)
    case ("d42")
        if(direction == "x") then
            call sbp_dx(Q42,lastnonzeroQ42,f,is,ie,js,je,n,is1,ie1,js1,je1,df)
        else !(direction == "y")
            call sbp_dy(Q42,lastnonzeroQ42,f,is,ie,js,je,n,is1,ie1,js1,je1,df)
        end if
    case default
        call parcomm_global%abort("sbp_mod, sbp_diff - unkonwn operator: "//opername)
    end select

end subroutine sbp_diff

subroutine sbp_dx(Q,last_nonzero,f,is,ie,js,je,nx,is1,ie1,js1,je1,dfdx)
    !dimensions of input function
    real(kind=8),     intent(in)    :: Q(:,:)
    integer(kind=4),  intent(in)    :: last_nonzero(:)
    integer(kind=4),  intent(in)    :: is, ie, js, je
    !global mesh dimension in x:
    integer(kind=4),  intent(in)    :: nx
    !output function dimensions:
    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1

    real(kind=8),     intent(in)    :: f(is:ie,js:je)
    real(kind=8),     intent(inout) :: dfdx(is1:ie1,js1:je1)

    integer(kind=4) i, j, mx, bs

    bs = size(Q,2)

    do j=js1,je1
        do i = max(is1,1),min(bs,ie1)
            mx = last_nonzero(i)
            dfdx(i,j) = sum(Q(1:mx,i)*f(1:mx,j))
        end do
        do i = max(bs+1,is1), min(ie1,nx-bs)
            dfdx(i,j) = (f(i-2,j)-8._8*f(i-1,j)+8._8*f(i+1,j)-f(i+2,j))/12.0_8
        end do
        do i=max(is1,nx-bs+1),min(ie1,nx)
            mx = last_nonzero(nx-i+1)
            dfdx(i,j) =-sum(Q(mx:1:-1,nx-i+1)*f(nx-mx+1:nx,j))
        end do
    end do
end subroutine sbp_dx

subroutine sbp_dy(Q,last_nonzero,f,is,ie,js,je,ny,is1,ie1,js1,je1,dfdy)
    !dimensions of input function
    real(kind=8),     intent(in)    :: Q(:,:)
    integer(kind=4),  intent(in)    :: last_nonzero(:)
    integer(kind=4),  intent(in)    :: is, ie, js, je
    !global mesh dimension in x:
    integer(kind=4),  intent(in)    :: ny
    !output function dimensions:
    integer(kind=4),  intent(in)    :: is1,ie1,js1,je1

    real(kind=8),     intent(in)    :: f(is:ie,js:je)
    real(kind=8),     intent(inout) :: dfdy(is1:ie1,js1:je1)

    integer(kind=4) i, j, my, bs

    bs = size(Q,2)

    do j = max(js1,1),min(bs,je1)
        do i=is1,ie1
            my = last_nonzero(j)
            dfdy(i,j) = sum(Q(1:my,j)*f(i,1:my))
        end do
    end do
    do j = max(bs+1,js1), min(je1,ny-bs)
        do i=is1,ie1
            dfdy(i,j) = (f(i,j-2)-8._8*f(i,j-1)+8._8*f(i,j+1)-f(i,j+2))/12.0_8
        end do
    end do
    do j=max(js1,ny-bs+1),min(je1,ny)
        do i=is1,ie1
            my = last_nonzero(ny-j+1)
            dfdy(i,j) =-sum(Q(my:1:-1,ny-j+1)*f(i,ny-my+1:ny))
        end do
    end do
end subroutine sbp_dy

end module sbp_mod
