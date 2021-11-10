module vertical_transform_mod

    use abstract_vertical_transform_mod, only : vertical_transform_t

    implicit none

    !basic case: eta = (z-h_surf)/(h_top-h_surf)
    !            z = h_surf+eta*(h_top-h_surf)
    type, public, extends(vertical_transform_t) :: vertical_transform_default_t
    contains
        procedure :: calc_z
        procedure :: calc_dz_deta
    end type

contains

    function calc_z(this, h_surf, h_top, eta) result(z)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: h_surf, h_top, eta
        real(kind=8) :: z

        z = h_surf+(h_top-h_surf)*eta
    end function calc_z

    function calc_dz_deta(this, h_surf, h_top, eta) result(dz_deta)
        class(vertical_transform_default_t), intent(in) :: this
        real(kind=8), intent(in) :: h_surf, h_top, eta
        real(kind=8) :: dz_deta

        dz_deta = (h_top-h_surf)
    end function calc_dz_deta

end module vertical_transform_mod
