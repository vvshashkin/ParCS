module halo_factory_mod

implicit none

contains

subroutine create_halo_procedure(halo,domain,halo_width,halo_type)
    use halo_mod,   only : halo_t
    use domain_mod, only : domain_t

    class(halo_t), allocatable, intent(out) :: halo
    class(domain_t),            intent(in)  :: domain
    integer(kind=4),            intent(in)  :: halo_width
    character(len=*),           intent(in)  :: halo_type

    if(halo_type=="A_default") then
        call create_A_default_halo_procedure(halo,domain,halo_width)
    else
        call domain%parcomm%abort("unknown halo_type in create_halo_procedure: "// &
                                   halo_type)
    end if
end

subroutine create_A_default_halo_procedure(halo,domain,halo_width)
    use halo_mod,           only : halo_t
    use domain_mod,         only : domain_t
    use halo_A_default_mod, only : halo_A_default_t
    use exchange_factory_mod,   only : create_symm_halo_exchange_A

    class(halo_t), allocatable, intent(out) :: halo
    class(domain_t),            intent(in)  :: domain
    integer(kind=4),            intent(in)  :: halo_width

    allocate(halo_A_default_t :: halo)
    select type(halo)
    type is (halo_A_default_t)
        halo%exch_halo = create_symm_halo_exchange_A( &
                   domain%partition, domain%parcomm, domain%topology, halo_width, 'full')
    end select
end

end module halo_factory_mod
