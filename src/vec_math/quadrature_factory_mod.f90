module quadrature_factory_mod

use abstract_quadrature_mod, only : quadrature_t
use default_quadrature_mod,  only : default_quadrature_t
use parcomm_mod,             only : parcomm_global
use domain_mod,              only : domain_t

implicit none

contains

subroutine create_quadrature(quadrature, quadrature_name, domain)

    class(quadrature_t), allocatable, intent(out) :: quadrature
    character(len=*),                 intent(in)  :: quadrature_name
    type(domain_t),                   intent(in)  :: domain

    if(quadrature_name == "default_quadrature") then
        quadrature = default_quadrature_t()
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_quadrature, quadrature_factory_mod")
    end if

end subroutine create_quadrature

end module quadrature_factory_mod
