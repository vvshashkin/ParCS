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

integer(kind=4), parameter :: ex(3,6) = reshape( (/  (/ 1, 0,  0/), &
                                                     (/ 0, 0,  1/), &
                                                     (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 1, 0,  0/)  /),  (/3 ,6/) )

integer(kind=4), parameter :: ey(3,6) = reshape( (/  (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 1, 0/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0, 1/)    /),  (/3 ,6/) )

integer(kind=4), parameter :: n(3,6) = reshape( (/   (/ 0, 0,  1/), &
                                                     (/-1, 0,  0/), &
                                                     (/ 0, 0, -1/), &
                                                     (/ 1, 0,  0/), &
                                                     (/ 0, 1,  0/), &
                                                     (/ 0,-1,  0/)    /),  (/3 ,6/) )

integer(kind=4), parameter ::  r(3,6) = reshape( (/  (/0, 0, 0/), &
                                                     (/1, 0, 0/), &
                                                     (/1, 0, 1/), &
                                                     (/0, 0, 1/), &
                                                     (/0, 0, 0/), &
                                                     (/0, 1, 0/)    /),  (/3 ,6/) )

integer(kind=4), protected, dimension(3,3,6) :: T, T2, T_right
integer(kind=4), protected :: Tinv(3,3,6)

contains

subroutine init_transform_matrices(Npoints)

    integer(kind=4), intent(in) :: Npoints

    integer(kind=4) :: panel_num

    do panel_num = 1, 6

        T2(1:3,1,panel_num) = ex(1:3, panel_num)
        T2(1:3,2,panel_num) = ey(1:3, panel_num)
        T2(1:3,3,panel_num) = r (1:3, panel_num)

        T_right(1:3,1,panel_num) = n (1:3, panel_num)
        T_right(1:3,2,panel_num) = ey(1:3, panel_num)
        T_right(1:3,3,panel_num) = r (1:3, panel_num)+ex(1:3, panel_num)

        T(1:3,1,panel_num) = ex(1:3, panel_num)
        T(1:3,2,panel_num) = ey(1:3, panel_num)
        T(1:3,3,panel_num) = (Npoints-1)*r (1:3, panel_num)+1 - ex(1:3, panel_num) - ey(1:3, panel_num)

        ! Tinv(1:3,1:3,panel_num) = matinv3(T(1:3, 1:3, panel_num))

    end do

    print*, T2(:,:,2)
end subroutine init_transform_matrices

subroutine apply_transform(transform, i, j, panel_num, x, y, z)

    integer(kind=4), intent(in)  :: transform(3,3,6)
    integer(kind=4), intent(in)  :: i, j, panel_num
    integer(kind=4), intent(out) :: x, y, z

    x = transform(1,1,panel_num)*i + transform(1,2,panel_num)*j + transform(1,3,panel_num)
    y = transform(2,1,panel_num)*i + transform(2,2,panel_num)*j + transform(2,3,panel_num)
    z = transform(3,1,panel_num)*i + transform(3,2,panel_num)*j + transform(3,3,panel_num)

end subroutine apply_transform

subroutine apply_T2_transform(i, j, panel_num, x, y, z)

    integer(kind=4), intent(in)  :: i, j, panel_num
    integer(kind=4), intent(out) :: x, y, z

    x = T2(1,1,panel_num)*i + T2(1,2,panel_num)*j + T2(1,3,panel_num)
    y = T2(2,1,panel_num)*i + T2(2,2,panel_num)*j + T2(2,3,panel_num)
    z = T2(3,1,panel_num)*i + T2(3,2,panel_num)*j + T2(3,3,panel_num)

end subroutine apply_T2_transform

subroutine apply_T_transform(i, j, panel_num, x, y, z)

    integer(kind=4), intent(in)  :: i, j, panel_num
    integer(kind=4), intent(out) :: x, y, z

    x = T(1,1,panel_num)*i + T(1,2,panel_num)*j + T(1,3,panel_num)
    y = T(2,1,panel_num)*i + T(2,2,panel_num)*j + T(2,3,panel_num)
    z = T(3,1,panel_num)*i + T(3,2,panel_num)*j + T(3,3,panel_num)

end subroutine apply_T_transform

subroutine change_basis(ex1, ey1, a1, i1, j1, ex2, ey2, a2, i2, j2)

    integer(kind=4), dimension(1:3), intent(in) :: ex1, ey1, a1
    integer(kind=4), dimension(1:3), intent(in) :: ex2, ey2, a2
    integer(kind=4), intent(in)  :: i1, j1
    integer(kind=4), intent(out) :: i2, j2



    i2 = dot_product((i1-1)*ex1 + (j1-1)*ey1 + a1 - a2, ex2) + 1
    j2 = dot_product((i1-1)*ex1 + (j1-1)*ey1 + a1 - a2, ey2) + 1


end subroutine change_basis

pure function matinv3(A) result(B)
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real(kind=4), intent(in) :: A(3,3)   !! Matrix
  real(kind=4)             :: B(3,3)   !! Inverse matrix

  real(kind=4)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

!
! subroutine transform_to_index(x,y,z, ix, iy, iz)
!
!
!
! end subroutine
!
! end subroutine transform_to_index

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

end module topology_mod
