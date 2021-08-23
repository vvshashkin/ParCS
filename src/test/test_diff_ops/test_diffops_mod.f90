module test_diffops_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use parcomm_mod,            only : parcomm_global

implicit none

type string_t
    character(:), allocatable :: str
end type string_t

type err_container_t
    type(string_t), allocatable :: keys(:)
    real(kind=8),   allocatable :: values(:)
end type err_container_t

private
public :: err_container_t, test_div, test_grad, test_laplace_spectre

real(kind=8), parameter :: some_const = 12.34567_8

contains

type(err_container_t) function test_div(N,div_oper_name,staggering) result(errs)

    use test_fields_mod,  only : set_vector_test_field, set_scalar_test_field, &
                                 solid_rot=>solid_rotation_field_generator, &
                                 cross_polar=>cross_polar_flow_generator, &
                                 cross_polar_div=>cross_polar_flow_div_generator
    use div_factory_mod,  only : create_div_operator
    use abstract_div_mod, only : div_operator_t

    integer(kind=4),  intent(in) :: N
    character(len=*), intent(in) :: div_oper_name, staggering
    !locals:
    integer(kind=4), parameter  :: nz = 3
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: u, v, div, div2, div_true
    type(domain_t)              :: domain
    class(div_operator_t), allocatable :: div_op

    call create_domain(domain, "cube", staggering, N, nz)
    call create_grid_field(u, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(v, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(div, ex_halo_width, 0, domain%mesh_p)
    call create_grid_field(div2, ex_halo_width, 0, domain%mesh_p)
    call create_grid_field(div_true, 0, 0, domain%mesh_p)

    call set_vector_test_field(u,v,solid_rot, domain%mesh_u, domain%mesh_v, &
                               0, "contravariant")

    div_op = create_div_operator(domain, div_oper_name)

    call div_op%calc_div(div, u,v,domain)
    call div_op%calc_div(div2,u,v,domain,some_const)

    call div2%update(-some_const,div,domain%mesh_p)
    if(div2%maxabs(domain%mesh_p,domain%parcomm)>1e-16) then
        call parcomm_global%abort("div_a2 test, wrong multiplier interface. Test failed!")
    end if

    allocate(errs%keys(4), errs%values(4))
    errs%keys(1)%str = "solid rotation linf"
    errs%keys(2)%str = "solid rotation l2"
    errs%keys(3)%str = "cross polar linf"
    errs%keys(4)%str = "cross polar l2"

    errs%values(1) = div%maxabs(domain%mesh_p,domain%parcomm)
    errs%values(2) = div%algebraic_norm2(domain%mesh_p,domain%parcomm)/real(N,8)

    call set_vector_test_field(u,v,cross_polar, domain%mesh_u, domain%mesh_v, &
                               0, "contravariant")
    call set_scalar_test_field(div_true,cross_polar_div, domain%mesh_p,0)
    call div_op%calc_div(div, u,v,domain)
    call div%update(-1.0_8,div_true,domain%mesh_p)

    errs%values(3) = div%maxabs(domain%mesh_p,domain%parcomm)
    errs%values(4) = div%algebraic_norm2(domain%mesh_p,domain%parcomm)/real(nz*N,8)
    !call stats(div,domain%mesh_p)

end function test_div

type(err_container_t) function test_grad(N,grad_oper_name,staggering) result(errs)

    use test_fields_mod,    only : set_vector_test_field, set_scalar_test_field, &
                                   xyz_f => xyz_scalar_field_generator, &
                                   xyz_grad => xyz_grad_generator
    use grad_factory_mod,   only : create_grad_operator
    use abstract_grad_mod,  only : grad_operator_t

    integer(kind=4),  intent(in) :: N
    character(len=*), intent(in) :: grad_oper_name, staggering
    !locals:
    integer(kind=4), parameter  :: nz = 3
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: gx, gy, gx1, gy1, f
    type(grid_field_t)          :: gx_true, gy_true
    type(domain_t)              :: domain
    class(grad_operator_t), allocatable :: grad_op

    call create_domain(domain, "cube", staggering, N, nz)
    call create_grid_field(gx, 0, 0, domain%mesh_u)
    call create_grid_field(gy, 0, 0, domain%mesh_v)
    call create_grid_field(gx1,0, 0, domain%mesh_u)
    call create_grid_field(gy1,0, 0, domain%mesh_v)
    call create_grid_field(gx_true, 0, 0, domain%mesh_u)
    call create_grid_field(gy_true, 0, 0, domain%mesh_v)
    call create_grid_field(f, ex_halo_width, 0, domain%mesh_p)

    call set_scalar_test_field(f,xyz_f, domain%mesh_p,0)
    call set_vector_test_field(gx_true, gy_true, xyz_grad, domain%mesh_u, domain%mesh_v, 0, "contravariant")

    grad_op = create_grad_operator(domain, grad_oper_name)
    call grad_op%calc_grad(gx,gy,f,domain)
    call grad_op%calc_grad(gx1,gy1,f,domain,some_const)

    call gx1%update(-some_const,gx,domain%mesh_u)
    call gy1%update(-some_const,gy,domain%mesh_v)

    if(gx1%maxabs(domain%mesh_u,domain%parcomm)+gy1%maxabs(domain%mesh_v,domain%parcomm)>1e-16) then
        call parcomm_global%abort("grad_a2 test, wrong multiplier interface. Test failed!")
    end if

    allocate(errs%keys(2), errs%values(2))
    errs%keys(1)%str = "xyz linf"
    errs%keys(2)%str = "xyz l2"

    call gx%update(-1.0_8,gx_true,domain%mesh_u)
    call gy%update(-1.0_8,gy_true,domain%mesh_v)
    errs%values(1) = gx%maxabs(domain%mesh_u,domain%parcomm)+ &
                     gy%maxabs(domain%mesh_v,domain%parcomm)
    errs%values(2) = gx%algebraic_norm2(domain%mesh_u,domain%parcomm)/real(nz*N,8)+&
                     gy%algebraic_norm2(domain%mesh_v,domain%parcomm)/real(nz*N,8)

    !call stats(gx,domain%mesh_u)

end function test_grad

subroutine test_laplace_spectre(div_operator_name, grad_operator_name, staggering)
    use test_fields_mod,   only : set_vector_test_field, solid_rot=>solid_rotation_field_generator
    use div_factory_mod,   only : create_div_operator
    use abstract_div_mod,  only : div_operator_t
    use grad_factory_mod,  only : create_grad_operator
    use abstract_grad_mod, only : grad_operator_t

    character(len=*), intent(in) :: div_operator_name, grad_operator_name, staggering

    integer(kind=4), parameter  :: nz = 1, nh = 12
    integer(kind=4), parameter  :: ex_halo_width = 8
    type(grid_field_t)          :: f, gx, gy, lap
    type(domain_t)              :: domain
    class(div_operator_t),  allocatable :: div_op
    class(grad_operator_t), allocatable :: grad_op
    real(kind=8), allocatable :: lapM(:,:)
    integer(kind=4) :: npts, ts, te, is, ie, js, je
    integer(kind=4) :: is1, ie1, js1, je1
    integer(kind=4) :: i, j, t, i1, j1, t1, ind, ind1

    if(parcomm_global%np>1) then
        call parcomm_global%print("ommiting laplace spectre test, works only in mpi n=1 mode")
        return
    end if

    call create_domain(domain, "cube", staggering, nh, nz)

    call create_grid_field(f,  ex_halo_width, 0, domain%mesh_p)
    call create_grid_field(gx, ex_halo_width, 0, domain%mesh_u)
    call create_grid_field(gy, ex_halo_width, 0, domain%mesh_v)
    call create_grid_field(lap, 0, 0, domain%mesh_p)

    grad_op = create_grad_operator(domain, grad_operator_name)
    div_op  = create_div_operator (domain, div_operator_name)

    ts = domain%mesh_p%ts
    te = domain%mesh_p%te
    npts = (te-ts+1)*nh**2
    allocate(lapM(npts,npts))

    call f%assign(0.0_8,domain%mesh_p)

    do t=ts, te
        is = domain%mesh_p%tile(t)%is
        ie = domain%mesh_p%tile(t)%ie
        js = domain%mesh_p%tile(t)%js
        je = domain%mesh_p%tile(t)%je
        do j=js,je
            do i=is,ie
                f%tile(t)%p(i,j,1) = 1.0_8
                call grad_op%calc_grad(gx,gy,f,domain)
                call div_op%calc_div(lap,gx,gy,domain)

                ind = (t-ts)*nh**2+(j-js)*nh+(i-is+1)

                do t1 = ts, te
                    is1 = domain%mesh_p%tile(t1)%is
                    ie1 = domain%mesh_p%tile(t1)%ie
                    js1 = domain%mesh_p%tile(t1)%js
                    je1 = domain%mesh_p%tile(t1)%je
                    do j1 = js1, je1
                        do i1= is1, ie1
                            ind1 = (t1-ts)*nh**2+(j1-js1)*nh+(i1-is1+1)
                            lapM(ind1,ind) = lap%tile(t1)%p(i1,j1,1)
                        end do
                    end do
                end do
                f%tile(t)%p(i,j,1) = 0.0_8
            end do
        end do
    end do

    call eigvals(lapM,npts)
end subroutine test_laplace_spectre

subroutine eigvals(H,M)
    integer(kind=4) M
    real(8) H(M,M)
    real(8) P(M,M)
    real(8) WR(M), WI(M)
    real(8) FV(M)
    complex(8) PP(M,M), P1(M,M), PR(M,M)
    complex(8) eval(M)
    integer IA(M)
    integer ierr, i, j
    complex(8), parameter :: im = (0._8,1._8)
    INTEGER IPVT(M)
    real(8) DET(2)
    complex(8) ZWK(2*M), ZDET(2)
    real(8) Rcond
    real(8), parameter :: zeps = 1e-14_8
    real(8) z

    call rg(M,M,H,WR,WI,1,P,IA,FV,ierr)

    print *, "Eigvals"
    print *, "max real", maxval(wr)
    print *, "min real", minval(wr)
    print *, "max imag", maxval(wi)
    print *, "min imag", minval(wi)

end subroutine eigvals

subroutine stats(f,mesh)
    use mesh_mod, only : mesh_t

    type(grid_field_t), intent(in) :: f
    type(mesh_t),       intent(in) :: mesh

    integer(kind=4) nx, ny, nz, t
    real(kind=8)  :: edge(4)
    real(kind=8)  :: corner(4)
    real(kind=8)  :: interior, interior_maxabs

    edge(1:4) = 0.0_8
    corner(1:4) = 0.0_8
    interior = 0.0_8
    interior_maxabs = 0.0_8

    do t = 1,6
        nx = mesh%tile(t)%ie
        ny = mesh%tile(t)%je
        nz = mesh%tile(t)%ke

        edge(1) = edge(1) + sum(abs(f%tile(t)%p(2:nx-1,1,:))) / (nx-2) / nz
        edge(2) = edge(2) + sum(abs(f%tile(t)%p(nx,2:ny-1,:))) / (ny-2) / nz
        edge(3) = edge(3) + sum(abs(f%tile(t)%p(2:nx-1,ny,:))) / (nx-2) / nz
        edge(4) = edge(4) + sum(abs(f%tile(t)%p(1,2:ny-1,:))) / (ny-2) / nz

        corner(1) = corner(1)+sum(abs(f%tile(t)%p(1,1,:)))/nz
        corner(2) = corner(2)+sum(abs(f%tile(t)%p(nx,1,:)))/nz
        corner(3) = corner(3)+sum(abs(f%tile(t)%p(nx,ny,:)))/nz
        corner(4) = corner(4)+sum(abs(f%tile(t)%p(1,ny,:)))/nz

        interior = interior + sum(abs(f%tile(t)%p(2:nx-1,2:ny-1,:))) / (nx-2)**2 / nz
        interior_maxabs = max(interior_maxabs, maxval(abs(f%tile(t)%p(2:nx-1,2:ny-1,:))))
    end do

    print '(A,4E15.7)', "edges: ", edge / 6.0_8
    print '(A,4E15.7)', "corner: ", corner / 6.0_8
    print *, "interior", interior/6.0, interior_maxabs
end subroutine stats

end module test_diffops_mod
