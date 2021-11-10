module abstract_vertical_transform_mod

    implicit none

    type, abstract :: vertical_transform_t
    contains
        procedure(calc_z_i), public, deferred :: calc_z
        procedure(calc_z_i), public, deferred :: calc_dz_deta
    end type

    abstract interface
        function calc_z_i(this, h_surf, h_top, eta) result(z)
            import vertical_transform_t
            class(vertical_transform_t), intent(in) :: this
            real(kind=8),                intent(in) :: h_surf, h_top, eta
            real(kind=8)                            :: z
        end function calc_z_i
    end interface

contains

end module abstract_vertical_transform_mod
