module parameters_NHlin_mod

use container_abstract_mod, only : model_parameters_abstract_t
use partition_mod,          only : partition_t
use mesh_mod,               only : mesh_t
use tile_mod,               only : tile_t

implicit none

type, extends(model_parameters_abstract_t) :: parameters_NHlin_t

    integer(kind=4)            :: nx     = 128
    integer(kind=4)            :: nz     = 10
    logical                    :: lcgrid = .true.

    real(kind=8)               :: HMAX   = 10e3_8
    real(kind=8)               :: dt     = 1._8
    real(kind=8)               :: radx   = 125._8 !Earth radius reduction

    integer(kind=4)            :: halo_width = 4

    !domain description
    integer(kind=4)            :: ts, te
    type(partition_t)          :: partition
    type(tile_t), allocatable  :: tiles(:)
    type(mesh_t), allocatable  :: mesh(:)

    !vertical grid
    real(kind=8)              :: dz
    real(kind=8), allocatable :: z(:), zh(:)

    !reference state
    real(kind=8), allocatable :: prex0(:), prex0dz(:)
    real(kind=8), allocatable :: theta0(:), theta0dz(:)

end type parameters_NHlin_t

contains

subroutine init_NHlin_parameters(params, namelist_str, myid, Np, master_id)

    use mesh_factory_mod,         only : create_equiangular_mesh

    class(parameters_NHlin_t), allocatable, intent(out) :: params
    character(:),              allocatable, intent(in)  :: namelist_str
    integer(kind=4),                        intent(in)  :: myid, Np, master_id


    integer(kind=4)            :: nx, nz
    logical                    :: lcgrid
    real(kind=8)               :: HMAX
    real(kind=8)               :: dt, radx

    namelist /dims/ nx, nz, lcgrid
    namelist /dyn/ HMAX, dt, radx

    integer(kind=4) ind

    allocate(parameters_NHlin_t :: params)

    if(allocated(namelist_str)) then
        !get defaults
        nx     = params%nx
        nz     = params%nz
        lcgrid = params%lcgrid
        HMAX   = params%HMAX
        dt     = params%dt
        radx   = params%radx
        !update with namelist values
        read(namelist_str, dims)
        read(namelist_str, dyn)
        !set namelist vals
        params%nx     = nx
        params%nz     = nz
        params%lcgrid = lcgrid
        params%HMAX   = HMAX
        params%dt     = dt
        params%radx   = radx
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

    call init_vertical_grid(params%nz,params%HMAX,namelist_str,params%z,params%zh,params%dz)
    call init_vertical_ref_profiles(params%nz, params%zh, params%z, params%prex0, params%prex0dz, &
                                                                    params%theta0, params%theta0dz)

    if(myid == master_id) then
        print *, "---Model parameters:"
        print *, "nx     =", params%nx
        print *, "nz     =", params%nz
        print *, "HMAX   =", params%HMAX
        print *, "dt     =", params%dt
        print *, "radx   =", params%radx
        print *, "lcgrid =", params%lcgrid
        print *, "--------------------"
        !print *, "zh", params%zh
        !print *, "z", params%z
        !print *, "prex", params%prex0
        !print *, "prexdz", params%prex0dz
        !print *, "theta", params%theta0
        !print *, "thetadz", params%theta0dz
    end if

end subroutine init_NHlin_parameters

subroutine init_vertical_grid(nz,ZTOP,namelist_str,z,zh,dz)

    integer(kind=4),             intent(in)  :: nz
    real(kind=8),                intent(in)  :: ZTOP
    character(:),   allocatable, intent(in)  :: namelist_str
    real(kind=8),   allocatable, intent(out) :: z(:), zh(:)
    real(kind=8),                intent(out) :: dz

    integer(kind=4) :: i

    dz = ZTOP / nz

    allocate(zh(0:nz))
    do i=0,nz
        zh(i) = i*dz
    end do

    allocate(z(1:nz))
    do i=1, nz
        z(i) = 0.5_8*(zh(i)+zh(i-1))
    end do

end subroutine init_vertical_grid

subroutine init_vertical_ref_profiles(nz, zh, z, prex0, prex0dz, theta0, theta0dz)

    use const_mod, only : grav, Cp

    integer(kind=4),           intent(in)  :: nz
    real(kind=8),              intent(in)  :: z(nz), zh(0:nz)
    real(kind=8), allocatable, intent(out) :: prex0(:), prex0dz(:)
    real(kind=8), allocatable, intent(out) :: theta0(:), theta0dz(:)

    real(kind=8), parameter :: RN      = 0.01_8
    real(kind=8), parameter :: RG      = grav**2/(RN**2*Cp)
    real(kind=8), parameter :: tconst  = 288._8

    integer k

    allocate(prex0(1:nz))
    allocate(prex0dz(1:nz))
    allocate(theta0(0:nz))
    allocate(theta0dz(0:nz))

    do k=1, nz
        prex0(k)    = ref_prex(z(k))
        prex0dz(k)  = ref_prexdz(z(k))
    end do

    do k=0, nz
        theta0(k)   = ref_theta(zh(k))
        theta0dz(k) = ref_thetadz(zh(k))
    end do

    contains

    real(8) function ref_prex(z) result(prex)
        real(8) z
        prex = (RG/tconst*exp(-RN**2/grav*z)+1._8-RG/tconst)
    end function ref_prex

    real(8) function ref_prexdz(z) result(prex)
        real(8) z
        prex =-RN**2*RG/(tconst*grav)*exp(-RN**2/grav*z)
    end function ref_prexdz

    real(8) function ref_theta(z) result(theta)
        real(8) z, temp
        theta = tconst*exp(RN**2/grav*z)
    end function ref_theta

    real(8) function ref_thetadz(z) result(theta)
        real(8) z, temp
        theta = RN**2/grav*tconst*exp(RN**2/grav*z)
    end function ref_thetadz

end subroutine init_vertical_ref_profiles

end module parameters_NHlin_mod
