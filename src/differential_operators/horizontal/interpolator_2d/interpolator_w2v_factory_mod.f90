module interpolator_w2v_factory_mod

use interpolator_w2v_mod, only : interpolator_w2v_t
use domain_mod,           only : domain_t
use sbp_factory_mod,      only : create_sbp_operator
use exchange_factory_mod, only : create_o_points_halo_exchange

implicit none

contains

subroutine create_w2v_interpolator(interpolator_w2v, sbp_w2v_interp_name, domain)

    type(interpolator_w2v_t), intent(out) :: interpolator_w2v
    character(len=*),         intent(in)  :: sbp_w2v_interp_name
    type(domain_t),           intent(in)  :: domain

    integer(kind=4) :: halo_width

    !WORKAROUND
    halo_width = 3

    interpolator_w2v%sbp_interp_w2v = create_sbp_operator(sbp_w2v_interp_name)

    interpolator_w2v%exchange = &
          create_o_points_halo_exchange(domain%partition, domain%parcomm, &
                                      domain%topology, halo_width, 'full')

end subroutine create_w2v_interpolator

end module interpolator_w2v_factory_mod
