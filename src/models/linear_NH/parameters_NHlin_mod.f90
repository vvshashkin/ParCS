module parameters_NHlin_mod

use container_abstract_mod, only : model_parameters_abstract_t
use partition_mod,          only : partition_t
use mesh_mod,               only : mesh_t
use tile_mod,               only : tile_t

implicit none

type, extends(model_parameters_abstract_t) :: parameters_NHlin_t

    integer(kind=4)            :: nx     = 128
    integer(kind=4)            :: nz     = 2
    logical                    :: lcgrid = .true.

    real(kind=8)               :: H0     = 1e3_8
    real(kind=8)               :: dt     = 300._8

    integer(kind=4)            :: halo_width = 8

    !domain description
    integer(kind=4)            :: ts, te
    type(partition_t)          :: partition
    type(tile_t), allocatable  :: tiles(:)
    type(mesh_t), allocatable  :: mesh(:)

end type parameters_NHlin_t

contains

subroutine init_NHlin_parameters(params, namelist_str, myid, Np, master_id)

    use mesh_factory_mod,         only : create_equiangular_mesh

    class(parameters_NHlin_t), allocatable, intent(out) :: params
    character(:),              allocatable, intent(in)  :: namelist_str
    integer(kind=4),                        intent(in)  :: myid, Np, master_id


    integer(kind=4)            :: nx, nz
    logical                    :: lcgrid
    real(kind=8)               :: H0
    real(kind=8)               :: dt

    namelist /dims/ nx, nz, lcgrid
    namelist /dyn/ H0, dt

    integer(kind=4) ind

    allocate(parameters_NHlin_t :: params)

    if(allocated(namelist_str)) then
        !get defaults
        nx     = params%nx
        nz     = params%nz
        lcgrid = params%lcgrid
        H0     = params%H0
        dt     = params%dt
        !update with namelist values
        read(namelist_str, dims)
        read(namelist_str, dyn)
        !set namelist vals
        params%nx     = nx
        params%nz     = nz
        params%lcgrid = lcgrid
        params%H0     = H0
        params%dt     = dt
    end if

    call params%partition%init(params%nx, params%nz, max(1,Np/6), myid, Np, strategy = 'default')

    params%ts = findloc(params%partition%proc_map, myid, dim=1)
    params%te = findloc(params%partition%proc_map, myid, back = .true., dim=1)

    allocate(params%tiles(params%ts:params%te))
    params%tiles(params%ts:params%te) = params%partition%tile(params%ts:params%te)

    allocate(params%mesh(params%ts:params%te))

    do ind = params%ts, params%te
        call create_equiangular_mesh(params%mesh(ind),                           &
                                     params%tiles(ind)%is, params%tiles(ind)%ie, &
                                     params%tiles(ind)%js, params%tiles(ind)%je, &
                                     params%tiles(ind)%ks, params%tiles(ind)%ke, &
                                     params%nx, params%halo_width,               &
                                     params%tiles(ind)%panel_number)
    end do

    if(myid == master_id) then
        print *, "---Model parameters:"
        print *, "nx=", params%nx
        print *, "nz=", params%nz
        print *, "H0=", params%H0
        print *, "dt=", params%dt
        print *, "lcgrid=", params%lcgrid
        print *, "--------------------"
    end if


end subroutine init_NHlin_parameters

end module parameters_NHlin_mod
