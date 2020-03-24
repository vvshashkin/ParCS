module NHlin_initial_cond_mod

implicit none

integer(kind=4)  :: test_case_num = 1
namelist /ini_cond/ test_case_num

private
public :: set_NHlin_initial_conditions

contains

subroutine set_NHlin_initial_conditions(stvec, namelist_str, params, &
                                        myid, master_id)

    use stvec_NHlin_mod,      only : stvec_NHlin_t
    use mesh_mod,             only : mesh_t
    use parameters_NHlin_mod, only : parameters_NHlin_t

    type(stvec_NHlin_t),       intent(inout) :: stvec
    character(:), allocatable, intent(in)    :: namelist_str
    type(parameters_NHlin_t),  intent(in)    :: params
    integer(kind=4),           intent(in)    :: myid, master_id

    integer(kind=4) namelist_read_stat

    if(allocated(namelist_str)) then
        read(namelist_str, ini_cond)
    end if

    if(myid == master_id) then
        print *, "---Initial conditions module initialization"
        print *, "test_case_num = ", test_case_num
        print *, "-------------------------------------------"
    end if

    if(test_case_num == 1) then
        call set_NHlin_gravity_wave(stvec, test_case_num, params%ts, params%te,  &
                                    params%mesh, params%radx, params%nz, params%zh, params%z)
    else
        call avost("NHlin model: unknown test case")
    end if

end subroutine set_NHlin_initial_conditions

subroutine set_NHlin_gravity_wave(stvec, test_case_num, ts, te, mesh, radx, nz, zh, z)

    use stvec_NHlin_mod, only : stvec_NHlin_t
    use mesh_mod,        only : mesh_t
    use const_mod,       only : pi, radz

    type(stvec_NHlin_t), intent(inout) :: stvec
    integer(kind=4),     intent(in)    :: test_case_num
    integer(kind=4),     intent(in)    :: ts, te
    type(mesh_t),        intent(in)    :: mesh(ts:te)
    real(kind=8),        intent(in)    :: radx
    integer(kind=4),     intent(in)    :: nz
    real(kind=8),        intent(in)    :: zh(0:nz), z(1:nz)

    real(kind=8), parameter :: x0 = 1._8/sqrt(2.0_8)
    real(kind=8), parameter :: y0 = 1._8/sqrt(2.0_8)
    real(kind=8), parameter :: z0 = 0._8
    real(kind=8), parameter :: r0 = 5000._8
    real(kind=8), parameter :: Lz = 20e3_8
    real(kind=8), parameter :: theta_max = 1.0_8

    integer(kind=4) ind
    integer(kind=4) i, j, k
    real(kind=8) dist

    do ind = ts, te
        stvec%prex(ind)%p(:,:,:) = 0._8
        stvec%u(ind)%p(:,:,:)    = 0._8
        stvec%v(ind)%p(:,:,:)    = 0._8
        stvec%w(ind)%p(:,:,:)    = 0._8

        do k = stvec%theta(ind)%ks, stvec%theta(ind)%ke
            do j = stvec%theta(ind)%js, stvec%theta(ind)%je
                do i = stvec%theta(ind)%is, stvec%theta(ind)%ie
                    dist = radz/radx*acos(mesh(ind)%rhx(i,j)*x0 + mesh(ind)%rhy(i,j)*y0 + &
                                          mesh(ind)%rhz(i,j)*z0)

                    stvec%theta(ind)%p(i,j,k) = theta_max*r0**2/(r0**2+dist**2) * &
                                                sin(2._8*pi*zh(k)/Lz)
                end do
            end do
        end do

    end do

end subroutine set_NHlin_gravity_wave

end module NHlin_initial_cond_mod
