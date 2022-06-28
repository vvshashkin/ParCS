module domain_factory_mod

use domain_mod,                only : domain_t
use topology_factory_mod,      only : create_topology
use cubed_sphere_topology_mod, only : cubed_sphere_topology_t
use metric_mod,                only : metric_t
use metric_factory_mod,        only : create_metric, create_metric_by_config
use orography_mod,             only : orography_t, orography_1mesh_t
use orography_factory_mod,     only : create_orography
use config_domain_mod,         only : config_domain_t
use mesh_factory_mod,          only : create_mesh
use parcomm_factory_mod,       only : create_parcomm
use parcomm_mod,               only : parcomm_global, parcomm_t
use mpi

implicit none

generic :: create_domain => create_domain_by_config, create_domain_by_arguments

contains

subroutine create_domain_by_arguments(domain, topology_type, staggering_type, nh, nz, &
                                      parcomm)

    use parcomm_mod,       only: parcomm_global, parcomm_t
    use config_domain_mod, only: config_domain_t

    type(domain_t),   intent(out) :: domain
    character(len=*), intent(in)  :: topology_type
    character(len=*), intent(in)  :: staggering_type
    integer(kind=4),  intent(in)  :: nh, nz
    type(parcomm_t),  optional, intent(in) :: parcomm

    type(config_domain_t) :: config_domain

    config_domain%N  = nh
    config_domain%Nz = nz
    config_domain%staggering_type     = staggering_type
    config_domain%vertical_staggering = "None"
    config_domain%metric_type         = "ecs"
    config_domain%topology_type       = topology_type
    config_domain%h_top               = 1.0_8
    call config_domain%config_metric%set_defaults()

    call create_domain_by_config(domain,config_domain,parcomm)

end subroutine create_domain_by_arguments

subroutine create_domain_by_config(domain, config, parcomm)

    type(domain_t),             intent(out) :: domain
    type(config_domain_t),      intent(in)  :: config
    type(parcomm_t),  optional, intent(in)  :: parcomm

    integer(kind=4) :: halo_width
    integer(kind=4) :: mpi_comm_local
    integer(kind=4) :: Nz_inter

    type(domain_t)  :: domain_2d

    domain%topology = create_topology(config%topology_type)

    call create_metric_by_config(domain%metric, domain%topology, &
                                 config%metric_type, config%config_metric)

    domain%horizontal_staggering = config%staggering_type

    mpi_comm_local = parcomm_global%comm_w
    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        !call create_parcomm(parcomm_glocal%comm_w, domain%parcomm)
        domain%parcomm = parcomm_global
    end if

    if(config%vertical_staggering == "CharneyPhilips") then
        Nz_inter = config%Nz+1
    else
        Nz_inter = config%Nz-1
    end if

    call domain%partition%init(config%N, config%Nz, Nz_inter, max(1,domain%parcomm%np/6), &
                               domain%parcomm%myid, domain%parcomm%Np,          &
                               config%staggering_type, strategy = 'default')

    if(config%is_orographic_curvilinear) then
        allocate(domain%domain_2d)
        call create_2d_subdomain(domain%domain_2d, config, domain%parcomm)
        call create_orography(domain%orography,config%orography_name,&
                              config%config_orography,domain%domain_2d,config%halo_width)
        call create_domain_meshes(domain,config,domain%orography)
    else
        call create_domain_meshes(domain,config)
    end if


end subroutine create_domain_by_config

subroutine create_2d_subdomain(domain, config, parcomm)

    use mesh_factory_mod,    only : create_mesh
    use parcomm_factory_mod, only : create_parcomm
    use parcomm_mod,         only : parcomm_global, parcomm_t

    type(domain_t),             intent(out) :: domain
    type(config_domain_t),      intent(in)  :: config
    type(parcomm_t),  optional, intent(in)  :: parcomm

    integer(kind=4) :: halo_width
    integer(kind=4) :: mpi_comm_local
    integer(kind=4) :: Nz_inter
    real(kind=8)    :: shift_zc, shift_zi, shift_xyz(3)

    type(domain_t)  :: domain_2d

    !have to be passed as an argument in future
    halo_width = 8

    domain%topology = create_topology(config%topology_type)

    call create_metric_by_config(domain%metric, domain%topology, &
                                 config%config_metric%metric_2d_type, config%config_metric)

    domain%horizontal_staggering = config%staggering_type

    if(present(parcomm)) then
        domain%parcomm = parcomm
    else
        domain%parcomm = parcomm_global
    end if

    call domain%partition%init(config%N, 1, 1, max(1,domain%parcomm%np/6), &
                               domain%parcomm%myid, domain%parcomm%Np,          &
                               config%staggering_type, strategy = 'default')

   call create_domain_meshes(domain,config)

end subroutine create_2d_subdomain

subroutine create_domain_meshes(domain,config, orography)

    type(domain_t),        intent(inout) :: domain
    type(config_domain_t), intent(in)    :: config
    type(orography_t),     intent(in), optional    :: orography

    real(kind=8)    :: shift_zc, shift_zi, shift_xyz(3)
    integer(kind=4) :: halo_width

    if(config%vertical_staggering == "CharneyPhilips") then
        shift_zc = 0.5_8; shift_zi = 0.0_8
    else
        shift_zc = 0.0_8; shift_zi = 0.5_8
    end if

    halo_width = config%halo_width

    shift_xyz(1:3) = [0.5_8, 0.5_8, shift_zc]
    call create_mesh(domain%mesh_o,  domain%partition, domain%metric, halo_width, &
                     config%h_top, domain%partition%tiles_o,  shift_xyz, orography%o)
    shift_xyz(1:3) = [0.0_8, 0.5_8, shift_zc]
    call create_mesh(domain%mesh_x, domain%partition, domain%metric, halo_width, &
                     config%h_top, domain%partition%tiles_x,  shift_xyz)
    shift_xyz(1:3) = [0.5_8, 0.0_8, shift_zc]
    call create_mesh(domain%mesh_y,  domain%partition, domain%metric, halo_width, &
                     config%h_top, domain%partition%tiles_y,  shift_xyz)
    shift_xyz(1:3) = [0.0_8, 0.0_8, shift_zc]
    call create_mesh(domain%mesh_xy, domain%partition, domain%metric, halo_width, &
                     config%h_top, domain%partition%tiles_xy, shift_xyz)
    if (config%vertical_staggering ==  "CharneyPhilips") then
        shift_xyz(1:3) = [0.5_8, 0.5_8, shift_zi]
        call create_mesh(domain%mesh_z,   domain%partition, domain%metric, halo_width, &
                         config%h_top, domain%partition%tiles_z, shift_xyz)
        shift_xyz(1:3) = [0.0_8, 0.0_8, shift_zi]
        call create_mesh(domain%mesh_xyz, domain%partition, domain%metric, halo_width, &
                         config%h_top, domain%partition%tiles_xyz, shift_xyz)
    end if

    select case(config%staggering_type)
    case ('A')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_o
        domain%mesh_v = domain%mesh_o
        domain%mesh_q = domain%mesh_o
    case ('Ah') !all degrees of freedom at corner points
        domain%mesh_p = domain%mesh_xy
        domain%mesh_u = domain%mesh_xy
        domain%mesh_v = domain%mesh_xy
        domain%mesh_q = domain%mesh_xy
    case ('C')
        domain%mesh_p = domain%mesh_o
        domain%mesh_u = domain%mesh_x
        domain%mesh_v = domain%mesh_y
        domain%mesh_q = domain%mesh_xy
    case ('Ch')
        domain%mesh_p = domain%mesh_xy
        domain%mesh_u = domain%mesh_y
        domain%mesh_v = domain%mesh_x
        domain%mesh_q = domain%mesh_o
    case default
        call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//config%staggering_type)
    end select
    select case(config%vertical_staggering)
    case("None")
        domain%mesh_w = domain%mesh_p
    case("CharneyPhilips")
            select case(config%staggering_type)
            case ('A','C')
                domain%mesh_w = domain%mesh_z
            case ('Ah','Ch')
                domain%mesh_w = domain%mesh_xyz
            case default
                call parcomm_global%abort("domain_factory_mod, unknown staggering type: "//config%staggering_type)
            end select
    case default
        call parcomm_global%abort("domain_factory_mod, unknown vertical staggering type: "//config%vertical_staggering)
    end select

end subroutine create_domain_meshes

end module domain_factory_mod
