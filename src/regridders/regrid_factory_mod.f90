module regrid_factory_mod

implicit none

contains

subroutine create_latlon_regrid(regrid_out, domain, Nlon, Nlat,interp_type, &
                                scalar_grid_type)

    use abstract_regrid_mod,    only : regrid_t
    use latlon_regrid_mod,      only : latlon_regrid_t, latlon_interp_t
    use grid_field_mod,         only : grid_field_t
    use domain_mod,             only : domain_t
    use parcomm_mod,            only : parcomm_global
    use halo_factory_mod,       only : create_halo_procedure
    use tile_mod,               only : tile_t

    class(regrid_t), allocatable, intent(out) :: regrid_out
    type(domain_t),               intent(in)  :: domain
    integer(kind=4),              intent(in)  :: Nlon, Nlat
    character(len=*),             intent(in)  :: interp_type
    character(len=*),             intent(in)  :: scalar_grid_type

    class(latlon_regrid_t), allocatable :: regrid
    type(tile_t) :: stencil_bounds
    integer(kind=4), parameter :: A_halo_width = 2

    allocate(regrid)

    regrid%Nlon = Nlon
    regrid%Nlat = Nlat
    regrid%scalar_grid_type = scalar_grid_type

    !Check:
    if(domain%partition%ts /= 1 .or. &
       domain%partition%te /= domain%partition%num_tiles*domain%partition%num_panels) then
       call parcomm_global%abort("cannot create latlon interpolator for distributed domains")
    end if

    if(regrid%scalar_grid_type == "Ah") then
        call stencil_bounds%init(is=1,ie=domain%mesh_o%tile(1)%nx+1,&
                                 js=1,je=domain%mesh_o%tile(1)%ny+1,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_scalar,domain,domain%mesh_xy, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    else if(regrid%scalar_grid_type == "A") then
        call create_halo_procedure(regrid%scalar_halo,domain,A_halo_width,"ECS_O")
        call stencil_bounds%init(is=1-A_halo_width,ie=domain%mesh_o%tile(1)%nx+A_halo_width,&
                                 js=1-A_halo_width,je=domain%mesh_o%tile(1)%ny+A_halo_width,&
                                 ks=domain%mesh_o%tile(1)%ks,ke=domain%mesh_o%tile(1)%ke)
        call create_latlon_interp(regrid%interp_scalar,domain,domain%mesh_o, &
                                  stencil_bounds, Nlat, Nlon, interp_type)
    end if

    call move_alloc(regrid,regrid_out)
end subroutine create_latlon_regrid

subroutine create_latlon_interp(interp,domain,mesh,stencil_bounds,Nlat,Nlon,interp_type)
    use latlon_regrid_mod,      only : latlon_interp_t
    use parcomm_mod,            only : parcomm_global
    use domain_mod,             only : domain_t
    use tile_mod,               only : tile_t
    use mesh_mod,               only : mesh_t

    use const_mod,              only : pi

    type(latlon_interp_t),  intent(out) :: interp
    type(domain_t),         intent(in)  :: domain
    type(mesh_t),           intent(in)  :: mesh
    type(tile_t),           intent(in)  :: stencil_bounds
    integer(kind=4),        intent(in)  :: Nlat, Nlon
    character(len=*),       intent(in)  :: interp_type

    real(kind=8)    :: lon, lat, hlon, hlat, r(3)
    real(kind=8)    :: alpha, beta
    integer(kind=4) :: panel_ind, i, j, ipanel, t, nx, ny
    integer(kind=4) :: iright
    real(kind=8)    :: hx, hy, alpha0, beta0, shift_alpha, shift_beta
    real(kind=8)    :: zdx, zdy

    interp%Nlon = Nlon
    interp%Nlat = Nlat
    interp%interp_type = interp_type
    interp%ks = mesh%tile(mesh%ts)%ks
    interp%ke = mesh%tile(mesh%ts)%ke

    select case(interp_type)
    case('linear')
        interp%stencil_width = 2
    case('cubic')
        interp%stencil_width = 4
    case default
        call parcomm_global%abort("regrid_factory_mod, create_latlon_interp"//&
                                  " - unknown interpolation type: "//interp_type)
    end select

    allocate(interp%wx(interp%stencil_width,Nlon,Nlat))
    allocate(interp%wy(interp%stencil_width,Nlon,Nlat))
    allocate(interp%tile_ind(Nlon,Nlat))
    allocate(interp%x_ind(Nlon,Nlat),interp%y_ind(Nlon,Nlat))

    hx = mesh%tile(mesh%ts)%hx
    hy = mesh%tile(mesh%ts)%hy
    alpha0 = mesh%tile(mesh%ts)%alpha_0
    beta0  = mesh%tile(mesh%ts)%beta_0
    shift_alpha = mesh%tile(mesh%ts)%shift_i
    shift_beta  = mesh%tile(mesh%ts)%shift_j
    nx = mesh%tile(mesh%ts)%nx
    ny = mesh%tile(mesh%ts)%ny

    hlon = 2._8*pi / real(Nlon,8)
    hlat = pi / (Nlat-1.0_8)

    do j = 1,Nlat
        lat = -0.5_8*pi+hlat*(j-1)
        do i = 1,Nlon
            lon = hlon*(i-1)
            r(1) = cos(lat)*cos(lon)
            r(2) = cos(lat)*sin(lon)
            r(3) = sin(lat)

            call domain%metric%transform_cartesian_to_native(panel_ind, alpha, beta, r)

            zdx = (alpha-alpha0)/hx+1-shift_alpha
            zdy = (beta -beta0) /hx+1-shift_beta
            !find leftmost point of stencil that lies in stencil bounds
            !we assume that (stencil_bounds%ie-stencil_bounds%is+1) > stencil_width
            !i.e. domain is wide enough
            iright = min(floor(zdx)+interp%stencil_width/2,stencil_bounds%ie)
            interp%x_ind(i,j) = max(iright-interp%stencil_width+1,stencil_bounds%is)
            iright = min(floor(zdy)+interp%stencil_width/2,stencil_bounds%je)
            interp%y_ind(i,j) = max(iright-interp%stencil_width+1,stencil_bounds%js)
            !save normalized displacement of point for future computation of interpolation weights
            interp%wx(1,i,j)  = zdx - interp%x_ind(i,j)
            interp%wy(1,i,j)  = zdy - interp%y_ind(i,j)
            !find tile
            do t = domain%partition%ts, domain%partition%te
                if(domain%partition%panel_map(t) == panel_ind   .and.&
                   mesh%tile(t)%is <= max(interp%x_ind(i,j),1)  .and.&
                   mesh%tile(t)%ie >= min(interp%x_ind(i,j),nx) .and.&
                   mesh%tile(t)%js <= max(interp%y_ind(i,j),1)  .and.&
                   mesh%tile(t)%je >= min(interp%y_ind(i,j),ny))  then
                   interp%tile_ind(i,j) = t
               end if
            end do
        end do
    end do

    if(interp_type == 'linear') then
        do j=1, Nlat
            do i=1, Nlon
                zdx = interp%wx(1,i,j)
                interp%wx(1,i,j) = 1._8-zdx
                interp%wx(2,i,j) = zdx
                zdy = interp%wy(1,i,j)
                interp%wy(1,i,j) = 1._8-zdy
                interp%wy(2,i,j) = zdy
            end do
        end do
    else if(interp_type == 'cubic') then
        do j=1, Nlat
            do i=1, Nlon
                zdx = interp%wx(1,i,j)
                zdy = interp%wy(1,i,j)
                interp%wx(1,i,j) =-(zdx-1._8)*(zdx-2._8)*(zdx-3._8) / 6._8
                interp%wx(2,i,j) = zdx*(zdx-2._8)*(zdx-3._8) / 2._8
                interp%wx(3,i,j) =-zdx*(zdx-1._8)*(zdx-3._8) / 2._8
                interp%wx(4,i,j) = zdx*(zdx-1._8)*(zdx-2._8) / 6._8
                interp%wy(1,i,j) =-(zdy-1._8)*(zdy-2._8)*(zdy-3._8) / 6._8
                interp%wy(2,i,j) = zdy*(zdy-2._8)*(zdy-3._8) / 2._8
                interp%wy(3,i,j) =-zdy*(zdy-1._8)*(zdy-3._8) / 2._8
                interp%wy(4,i,j) = zdy*(zdy-1._8)*(zdy-2._8) / 6._8
            end do
        end do
    end if
end

end module regrid_factory_mod
