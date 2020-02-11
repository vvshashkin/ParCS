module parameters_iomega_mod

use container_abstract_mod, only: model_parameters_abstract_t

implicit none

type, extends(model_parameters_abstract_t) :: parameters_iomega_t

end type parameters_iomega_t

end module parameters_iomega_mod
