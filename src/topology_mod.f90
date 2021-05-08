module topology_mod

implicit none

!!             ......
!!            |      :
!!            |  6   :
!!            |______:
!!   ......    ......   ......    ......
!!  |      :  |      : |      :  |      :
!!  |   4  :  |   1  : |   2  :  |   3  :
!!  |______:  |______: |______:  |______:
!!             ______
!!            |      :
!!            |  5   :
!!            |...,,,:
!!

integer(kind=4), parameter :: ex(3,6) = reshape( (/  (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 0, 0,  1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 1, 0,  0/)  /),  (/3 ,6/) )

integer(kind=4), parameter :: ey(3,6) = reshape( (/  (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0,-1/)    /),  (/3 ,6/) )

integer(kind=4), parameter :: n(3,6) = reshape( (/   (/ 0, 0,  1/), &
                                                     (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 0, 1,  0/), &
                                                     (/ 0,-1,  0/)    /),  (/3 ,6/) )

integer(kind=4), parameter ::  r(3,6) = reshape( (/  (/1, 0, 0/), &
                                                     (/1, 0, 1/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0, 0/), &
                                                     (/0, 0, 0/), &
                                                     (/0, 1, 1/)    /),  (/3 ,6/) )
contains

subroutine calc_xyz_cords(panel_ind, i, j, npoints, ix, iy, iz)

    integer(kind=4), intent(in)  :: panel_ind, i, j, npoints
    integer(kind=4), intent(out) :: ix, iy, iz

    ix = (i-1)*ex(1, panel_ind) + (j-1)*ey(1, panel_ind) + (npoints-1)*r(1,panel_ind) + 1
    iy = (i-1)*ex(2, panel_ind) + (j-1)*ey(2, panel_ind) + (npoints-1)*r(2,panel_ind) + 1
    iz = (i-1)*ex(3, panel_ind) + (j-1)*ey(3, panel_ind) + (npoints-1)*r(3,panel_ind) + 1

end subroutine calc_xyz_cords

subroutine transform_index(pn_out, pn_in, Npoints, i_in, j_in, i_out, j_out, i_step, j_step, first_dim_index)

    integer(kind=4),  intent(in)            :: pn_out, pn_in, Npoints, i_in, j_in
    integer(kind=4),  intent(out)           :: i_out , j_out
    integer(kind=4),  intent(out), optional :: i_step, j_step  !Determines descending or ascending order of loop
    character(len=1), intent(out), optional :: first_dim_index !Determines first loop index in horizontal direction

    integer(kind=4) :: ex_local(1:3), ey_local(1:3), r_local(1:3)
    integer(kind=4) :: ex_ex, ex_ey, ex_r, ey_ex, ey_ey, ey_r
    logical :: is_opposite, is_right, is_left, is_bottom, is_top


    is_right    = (dot_product(ex(1:3,pn_out), n(1:3,pn_in)) == -1)
    is_left     = (dot_product(ex(1:3,pn_out), n(1:3,pn_in)) ==  1)
    is_bottom   = (dot_product(ey(1:3,pn_out), n(1:3,pn_in)) ==  1)
    is_top      = (dot_product(ey(1:3,pn_out), n(1:3,pn_in)) == -1)
    is_opposite = (dot_product( n(1:3,pn_out), n(1:3,pn_in)) == -1)

    if (is_right) then

        ex_local = n (1:3, pn_out)
        ey_local = ey(1:3, pn_out)
        r_local  = r (1:3, pn_out) + ex(1:3, pn_out)

        ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
        ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
        ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

        ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
        ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
        ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

        i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (Npoints-1)*ex_r + 1 + Npoints
        j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (Npoints-1)*ey_r + 1

    else if (is_left) then

        ex_local = -n (1:3, pn_out)
        ey_local =  ey(1:3, pn_out)
        r_local  =  r (1:3, pn_out) - ex_local(1:3)

        ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
        ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
        ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

        ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
        ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
        ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

        i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (Npoints-1)*ex_r + 1 - Npoints
        j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (Npoints-1)*ey_r + 1

    else if (is_bottom) then

        ex_local =  ex(1:3, pn_out)
        ey_local = -n(1:3, pn_out)
        r_local  =  r (1:3, pn_out) + n(1:3, pn_out)

        ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
        ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
        ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

        ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
        ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
        ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

        i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (Npoints-1)*ex_r + 1
        j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (Npoints-1)*ey_r + 1 - Npoints

    else if (is_top) then

        ex_local =  ex(1:3, pn_out)
        ey_local =  n(1:3, pn_out)
        r_local  =  r(1:3, pn_out) +ey(1:3,pn_out)

        ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
        ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
        ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

        ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
        ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
        ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

        i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (Npoints-1)*ex_r + 1
        j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (Npoints-1)*ey_r + 1 + Npoints

    else if (is_opposite) then

        ex_local = -ex(1:3, pn_out)
        ey_local =  ey(1:3, pn_out)
        r_local  =  r (1:3, pn_out) + ex(1:3, pn_out) + n(1:3, pn_out)

        ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
        ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
        ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

        ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
        ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
        ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

        i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (Npoints-1)*ex_r + 1 + 2* Npoints
        j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (Npoints-1)*ey_r + 1

    end if

    if (present(i_step).and.present(j_step).and. present(first_dim_index) ) then
        if ( ex_ex /= 0 ) then
            first_dim_index = 'i'
            i_step = ex_ex
            j_step = ey_ey
        else
            first_dim_index = 'j'
            i_step = ex_ey
            j_step = ey_ex
        end if
    end if

end subroutine transform_index
subroutine transform_index_new(pn_in, i_in, j_in, nx, ny, pn_out, i_out, j_out, i_step, j_step, first_dim_index)

    integer(kind=4),  intent(in)            :: pn_out, pn_in, nx, ny, i_in, j_in
    integer(kind=4),  intent(out)           :: i_out , j_out
    integer(kind=4),  intent(out), optional :: i_step, j_step  !Determines descending or ascending order of loop
    character(len=1), intent(out), optional :: first_dim_index !Determines first loop index in horizontal direction

    integer(kind=4) :: ex_local(1:3), ey_local(1:3), r_local(1:3), shift(2)
    integer(kind=4) :: ex_ex, ex_ey, ex_r, ey_ex, ey_ey, ey_r
    logical :: is_opposite, is_right, is_left, is_bottom, is_top

    is_right    = (dot_product(ex(1:3,pn_out), n(1:3,pn_in)) == -1)
    is_left     = (dot_product(ex(1:3,pn_out), n(1:3,pn_in)) ==  1)
    is_bottom   = (dot_product(ey(1:3,pn_out), n(1:3,pn_in)) ==  1)
    is_top      = (dot_product(ey(1:3,pn_out), n(1:3,pn_in)) == -1)
    is_opposite = (dot_product( n(1:3,pn_out), n(1:3,pn_in)) == -1)

    if (is_right) then
        ex_local = n (1:3, pn_out)
        ey_local = ey(1:3, pn_out)
        r_local  = r (1:3, pn_out) + ex(1:3, pn_out)
        shift = [nx,0]
    else if (is_left) then
        ex_local = -n (1:3, pn_out)
        ey_local =  ey(1:3, pn_out)
        r_local  =  r (1:3, pn_out) - ex_local(1:3)
        shift = [-nx,0]
    else if (is_bottom) then
        ex_local =  ex(1:3, pn_out)
        ey_local = -n(1:3, pn_out)
        r_local  =  r (1:3, pn_out) + n(1:3, pn_out)
        shift = [0,-ny]
    else if (is_top) then
        ex_local =  ex(1:3, pn_out)
        ey_local =  n(1:3, pn_out)
        r_local  =  r(1:3, pn_out) +ey(1:3,pn_out)
        shift = [0, ny]
    else if (is_opposite) then
        ex_local = -ex(1:3, pn_out)
        ey_local =  ey(1:3, pn_out)
        r_local  =  r (1:3, pn_out) + ex(1:3, pn_out) + n(1:3, pn_out)
        shift = [2*nx,0]
    end if

    ex_ex = dot_product(ex(1:3, pn_in), ex_local(1:3))
    ex_ey = dot_product(ey(1:3, pn_in), ex_local(1:3))
    ex_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ex_local(1:3))

    ey_ex = dot_product(ex(1:3, pn_in), ey_local(1:3))
    ey_ey = dot_product(ey(1:3, pn_in), ey_local(1:3))
    ey_r  = dot_product(r(1:3, pn_in) - r_local(1:3), ey_local(1:3))

    i_out = (i_in-1)*ex_ex + (j_in-1)*ex_ey + (nx-1)*ex_r + 1 + shift(1)
    j_out = (i_in-1)*ey_ex + (j_in-1)*ey_ey + (ny-1)*ey_r + 1 + shift(2)

    if (present(i_step).and.present(j_step).and. present(first_dim_index) ) then
        if ( ex_ex /= 0 ) then
            first_dim_index = 'i'
            i_step = ex_ex
            j_step = ey_ey
        else
            first_dim_index = 'j'
            i_step = ex_ey
            j_step = ey_ex
        end if
    end if

end subroutine transform_index_new
end module topology_mod
