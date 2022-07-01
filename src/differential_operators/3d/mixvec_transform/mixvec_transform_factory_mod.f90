module mixvec_transform_factory_mod

use abstract_mixvec_transform_mod,   only : mixvec_transform_t
use mixvec_transform_colocated_mod,  only : mixvec_transform_colocated_t
use domain_mod,                      only : domain_t
use parcomm_mod,                     only : parcomm_global

implicit none

contains

subroutine create_mixvec_transform(mixvec_transform,mixvec_transform_name,domain)
    class(mixvec_transform_t), allocatable, intent(out) :: mixvec_transform
    character(len=*),                       intent(in)  :: mixvec_transform_name
    type(domain_t),                         intent(in)  :: domain

    select case(mixvec_transform_name)
    case("mixvec_colocated")
        mixvec_transform = mixvec_transform_colocated_t()
    case default
        call parcomm_global%abort("mixvec_transform_factory_mod, incorrect mixvec_transform_name: "//mixvec_transform_name)
    end select
end subroutine create_mixvec_transform

end module mixvec_transform_factory_mod
