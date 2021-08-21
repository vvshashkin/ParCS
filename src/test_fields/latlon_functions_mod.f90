module latlon_functions_mod
    implicit none
contains
    real(kind=8) pure function sin_phi(x,y,z) result(s)
        real(kind=8), intent(in) :: x,y,z
        s = z / sqrt(x**2+y**2+z**2)
    end function sin_phi

    real(kind=8) pure function cos_phi(x,y,z) result(c)
        real(kind=8), intent(in) :: x,y,z
        c = sqrt((x**2+y**2)/(x**2+y**2+z**2))
    end function cos_phi

    real(kind=8) pure function cos_lam(x,y,z) result(c)
        real(kind=8), intent(in) :: x,y,z
        c = x / max(sqrt(x**2+y**2),1e-14)
    end function cos_lam

    real(kind=8) pure function sin_lam(x,y,z) result(s)
        real(kind=8), intent(in) :: x,y,z
        !s = cos_lam(x,y,z)
        !s = sqrt(1.0_8-s**2)
        s = y / max(sqrt(x**2+y**2),1e-14)
        if(x**2+y**2 == 0) s = 1.0_8
    end function sin_lam
end module latlon_functions_mod
