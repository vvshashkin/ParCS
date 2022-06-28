module config_orography_mod

use config_mod, only : config_t

implicit none

type, extends(config_t) :: config_test_orography_t
    real(kind=8) :: h = 1.0_8
end type

end module config_orography_mod
