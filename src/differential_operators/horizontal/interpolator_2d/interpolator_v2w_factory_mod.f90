module interpolator_v2w_factory_mod

use interpolator_v2w_mod, only : interpolator_v2w_t
use domain_mod,           only : domain_t
use sbp_factory_mod,      only : create_sbp_operator
use exchange_factory_mod, only : create_symmetric_halo_vec_exchange_C

implicit none

contains

subroutine create_v2w_interpolator(interpolator_v2w, sbp_v2w_interp_name, domain)

    type(interpolator_v2w_t), intent(out) :: interpolator_v2w
    character(len=*),         intent(in)  :: sbp_v2w_interp_name
    type(domain_t),           intent(in)  :: domain

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 3

    interpolator_v2w%sbp_interp_v2w = create_sbp_operator(sbp_v2w_interp_name)

    interpolator_v2w%exchange = &
          create_symmetric_halo_vec_exchange_C(domain%partition, domain%parcomm, &
                                               domain%topology, halo_width, 'full')

end subroutine create_v2w_interpolator

end module interpolator_v2w_factory_mod
