module orography_factory_mod

use orography_mod, only : orography_t

implicit none

contains

subroutine create_orography(orography)
    type(orography_t), intent(out) :: orography

    print *, "Hello"
end subroutine

end module orography_factory_mod
