module test_halo_mod

implicit none

contains

subroutine test_halo()

    use domain_mod,         only : domain_t
    use domain_factory_mod, only : create_domain
    use halo_mod,           only : halo_t
    use halo_factory_mod,   only : create_halo_procedure
    use grid_field_mod,         only : grid_field_t
    use grid_field_factory_mod, only : create_grid_field

    type(domain_t) domain
    class(halo_t), allocatable :: halo

    integer(kind=4), parameter :: nh = 120,nz=2,halo_width=3

    type(grid_field_t) :: f
    integer(kind=4) :: t, i, j, k

    call create_domain(domain, "cube", 'A', nh, nz)
    call create_halo_procedure(halo,domain,halo_width,"A_default")

    call create_grid_field(f, halo_width, 0, domain%mesh_p)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    f%tile(t)%p(i,j,k) = (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
                end do
            end do
        end do
    end do

    call halo%get_halo_scalar(f,domain%parcomm,halo_width)

end subroutine test_halo

end module test_halo_mod
