module sbp_operator_mod

use grid_field_mod, only : tile_field_t
use tile_mod,       only : tile_t
use parcomm_mod,    only : parcomm_global

implicit none

!Implements application of 1D SBP operator (difference or interpolation)

type sbp_operator_t
    real(kind=8), allocatable    :: W_edge_l(:,:) !operator weights at left edge
    real(kind=8), allocatable    :: W_edge_r(:,:) !operator weights at right edge
    integer(kind=4), allocatable :: edge_last_l(:)  !last nonzero weight in W_edge_l
    integer(kind=4), allocatable :: edge_first_shift_r(:) !first nonzero weight in W_edge_r (minus n)
    real(kind=8), allocatable    :: W_in(:)         !inner stencil weights
    integer(kind=4)              :: in_shift        !shift of inner stencil first element with respect to
                                                    !the index of calculated element
    real(kind=8), allocatable    :: Al_out(:,:), Al_in(:,:) !Mass matrices at output and input grids left side
    real(kind=8), allocatable    :: Ar_out(:,:), Ar_in(:,:) !Mass matrices at output and input grids right side
    integer(kind=4)              :: dnx  ! difference of the dimensions of output and input grids
                                         ! in the direction of operator application
                                         ! nx_output - nx_input

    contains
    !specialization on different types of input/output arguments
    !do we really need it?
    !procedure, public :: apply_gf_to_gf_mesh    => apply_sbp_gf_to_gf_mesh
    !procedure, public :: apply_gf_to_gf_tile    => apply_sbp_gf_to_gf_tile

    !Grid field in -> array out
    procedure, public :: apply_gf_to_array      => apply_sbp_gf_to_array
    !array in -> array out
    procedure, public :: apply_array_to_array   => apply_sbp_array_to_array
    generic :: apply => apply_gf_to_array, apply_array_to_array

end type sbp_operator_t

contains

subroutine apply_sbp_gf_to_array(this, a_out, work_tile, out_tile, &
                                                nx_out_grid, direction, f_in)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    !bounding_tile = output arrat bounds
    type(tile_t),          intent(in) :: work_tile, out_tile
    integer(kind=4),       intent(in) :: nx_out_grid
    character(len=*),      intent(in) :: direction
    type(tile_field_t),    intent(in) :: f_in
    !output:
    real(kind=8), intent(out) :: a_out(out_tile%is:out_tile%ie, &
                                       out_tile%js:out_tile%je, &
                                       out_tile%ks:out_tile%ke)


    integer(kind=4) :: nx_in_grid, k, f_is, f_ie, f_js, f_je

    nx_in_grid = nx_out_grid-this%dnx

    do k=work_tile%ks, work_tile%ke
        f_is = f_in%is; f_ie = f_in%ie
        f_js = f_in%js; f_je = f_in%je

        call sbp_apply(a_out(out_tile%is:out_tile%ie, out_tile%js:out_tile%je, k),&
                       out_tile%is,  out_tile%ie,  out_tile%js,  out_tile%je,     &
                       work_tile%is, work_tile%ie, work_tile%js, work_tile%je,    &
                       f_in%p(f_is:f_ie,f_js:f_je,k), f_is, f_ie, f_js, f_je,     &
                       nx_in_grid, nx_out_grid,                                   &
                       this%W_edge_l, this%edge_last_l,                           &
                       this%W_edge_r, this%edge_first_shift_r, &
                       this%W_in, this%in_shift, direction)
    end do

end subroutine apply_sbp_gf_to_array

subroutine apply_sbp_array_to_array(this, a_out, work_tile, out_tile, &
                                      nx_out_grid, direction, a_in, in_tile)
    class(sbp_operator_t), intent(in) :: this
    !work_tile = indices where operator action will be calculated
    !bounding_tile = output arrat bounds
    type(tile_t),          intent(in) :: work_tile, out_tile
    integer(kind=4),       intent(in) :: nx_out_grid
    character(len=*),      intent(in) :: direction
    type(tile_t),          intent(in) :: in_tile
    real(kind=8),          intent(in) :: a_in(in_tile%is:in_tile%ie, &
                                              in_tile%js:in_tile%je, &
                                              in_tile%ks:in_tile%ke)
    !output:
    real(kind=8), intent(out) :: a_out(out_tile%is:out_tile%ie, &
                                       out_tile%js:out_tile%je, &
                                       out_tile%ks:out_tile%ke)


    integer(kind=4) :: nx_in_grid, k, f_is, f_ie, f_js, f_je

    nx_in_grid = nx_out_grid-this%dnx

    do k=work_tile%ks, work_tile%ke
        f_is = in_tile%is; f_ie = in_tile%ie
        f_js = in_tile%js; f_je = in_tile%je

        call sbp_apply(a_out(out_tile%is:out_tile%ie, out_tile%js:out_tile%je, k),&
                       out_tile%is,  out_tile%ie,  out_tile%js,  out_tile%je,     &
                       work_tile%is, work_tile%ie, work_tile%js, work_tile%je,    &
                       a_in(f_is:f_ie,f_js:f_je,k), f_is, f_ie, f_js, f_je,     &
                       nx_in_grid, nx_out_grid,                                   &
                       this%W_edge_l, this%edge_last_l,                           &
                       this%W_edge_r, this%edge_first_shift_r, &
                       this%W_in, this%in_shift, direction)
    end do

end subroutine apply_sbp_array_to_array

subroutine sbp_apply(df,is1,ie1,js1,je1,is,ie,js,je,                          &
                     f,is0,ie0,js0,je0,nsrc,ntarg,                            &
                     Ql,l_last_nonzero,Qr,r_first_nonzero_shift,Qin,inshift,direction)
    !dimensions of input function
    real(kind=8),     intent(in)    :: Ql(:,:)
    integer(kind=4),  intent(in)    :: l_last_nonzero(:)
    real(kind=8),     intent(in)    :: Qr(:,:)
    integer(kind=4),  intent(in)    :: r_first_nonzero_shift(:)
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

    real(kind=8),     intent(inout) :: df(is1:ie1,js1:je1)

    integer(kind=4) :: i, ir, j, jr
    integer(kind=4) :: mx, my, nwst, nwr
    integer(kind=4) :: l_edge_size, r_edge_size, inwidth, inend

    l_edge_size = size(Ql,2)
    r_edge_size = size(Qr,2)
    nwr = size(Qr,1)
    inwidth = size(Qin)
    inend = inwidth+inshift-1

    select case(direction)
    case('x')

    do j=js,je
        do i = max(is,1),min(l_edge_size,ie)
            mx = l_last_nonzero(i)
            df(i,j) = sum(Ql(1:mx,i)*f(1:mx,j))
        end do
        do i = max(l_edge_size+1,is), min(ie,ntarg-r_edge_size)
            df(i,j) = sum(Qin(1:inwidth)*f(i+inshift:i+inend,j))
        end do
        do i=max(is,ntarg-r_edge_size+1),min(ie,ntarg)
            ir = i - (ntarg-r_edge_size)
            mx = nsrc+r_first_nonzero_shift(ir)
            nwst = nwr+r_first_nonzero_shift(ir)
            df(i,j) = sum(Qr(nwst:nwr,ir)*f(mx:nsrc,j))
        end do
    end do

    case('y')

    do j = max(js,1),min(l_edge_size,je)
        do i=is,ie
            my = l_last_nonzero(j)
            df(i,j) = sum(Ql(1:my,j)*f(i,1:my))
        end do
    end do
    do j = max(l_edge_size+1,js), min(je,ntarg-r_edge_size)
        do i=is,ie
            df(i,j) = sum(Qin(1:inwidth)*f(i,j+inshift:j+inend))
        end do
    end do
    do j=max(js,ntarg-r_edge_size+1),min(je,ntarg)
        do i=is,ie
            jr = j - (ntarg-r_edge_size)
            my = nsrc+r_first_nonzero_shift(jr)
            nwst = nwr+r_first_nonzero_shift(jr)
            df(i,j) = sum(Qr(nwst:nwr,jr)*f(i,my:nsrc))
        end do
    end do

    case default
    call parcomm_global%abort("unknown direction in sbp_apply: "//direction)
    end select
end subroutine sbp_apply

end module sbp_operator_mod
