module swlin_initial_cond_mod

implicit none

integer(kind=4)  :: test_case_num = 1
namelist /ini_cond/ test_case_num

private
public :: set_swlin_initial_conditions

contains

subroutine set_swlin_initial_conditions(stvec, namelist_str, ts, te, mesh, &
                                        myid, master_id)

    use stvec_swlin_mod, only : stvec_swlin_t
    use mesh_mod,        only : mesh_t

    type(stvec_swlin_t),       intent(inout) :: stvec
    character(:), allocatable, intent(in)    :: namelist_str
    integer(kind=4),           intent(in)    :: ts, te
    type(mesh_t),              intent(in)    :: mesh(ts:te)
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
        call set_swlin_gauss_hill(stvec, test_case_num, ts, te, mesh)
    else
        call avost("swlin model: unknown test case")
    end if

end subroutine set_swlin_initial_conditions

subroutine set_swlin_gauss_hill(stvec, test_case_num, ts, te, mesh)

    use stvec_swlin_mod, only : stvec_swlin_t
    use mesh_mod,        only : mesh_t

    type(stvec_swlin_t), intent(inout) :: stvec
    integer(kind=4),     intent(in)    :: test_case_num
    integer(kind=4),     intent(in)    :: ts, te
    type(mesh_t),        intent(in)    :: mesh(ts:te)

    real(kind=8), parameter :: x0 = 0._8, y0 = 0._8, z0 = 1._8
    real(kind=8), parameter :: r0 = 0.25_8
    real(kind=8), parameter :: hmax = 1.0_8

    integer(kind=4) ind
    integer(kind=4) i, j
    real(kind=8) dist

    do ind = ts, te
        stvec%u(ind)%p(:,:,:) = 0._8
        stvec%v(ind)%p(:,:,:) = 0._8

        do j = stvec%h(ind)%js, stvec%h(ind)%je
            do i = stvec%h(ind)%is, stvec%h(ind)%ie
                dist = acos(mesh(ind)%rhx(i,j)*x0+mesh(ind)%rhy(i,j)*y0+ &
                                                  mesh(ind)%rhz(i,j)*z0)

                stvec%h(ind)%p(i,j,1) = hmax*exp(-(dist/r0)**2)
            end do
        end do

    end do

end subroutine set_swlin_gauss_hill

end module swlin_initial_cond_mod
