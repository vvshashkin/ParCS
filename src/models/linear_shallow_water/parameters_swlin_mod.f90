module parameters_swlin_mod

use container_abstract_mod, only: model_parameters_abstract_t

implicit none

type, extends(model_parameters_abstract_t) :: parameters_swlin_t

end type parameters_swlin_t

end module parameters_swlin_mod
