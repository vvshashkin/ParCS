module interpolator_h2v_factory_mod

use interpolator_h2v_mod, only : interpolator_h2v_t
use domain_mod,           only : domain_t
use sbp_factory_mod,      only : create_sbp_operator
use exchange_factory_mod, only : create_symm_halo_exchange_A

implicit none

contains

subroutine create_h2v_interpolator(interpolator_h2v, sbp_h2v_interp_name, domain)

    type(interpolator_h2v_t), intent(out) :: interpolator_h2v
    character(len=*),         intent(in)  :: sbp_h2v_interp_name
    type(domain_t),           intent(in)  :: domain

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 3

    interpolator_h2v%sbp_interp_h2v = create_sbp_operator(sbp_h2v_interp_name)

    interpolator_h2v%exchange = &
          create_symm_halo_exchange_A(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')

end subroutine create_h2v_interpolator

end module interpolator_h2v_factory_mod
