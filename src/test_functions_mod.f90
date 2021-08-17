module test_fields_mod
use grid_field_mod, only : grid_field_t
use mesh_mod,       only : mesh_t

implicit none

private
public :: set_scalar_test_field, set_vector_test_field
public :: xyz_scalar_field_generator_t, xyz_scalar_field_generator
public :: solid_rotation_field_generator_t, solid_rotation_field_generator
public :: xyz_grad_generator_t, xyz_grad_generator

!!!!!!!!!!!!!Abstract scalar and vector fields generators
type, public, abstract :: scalar_field_generator_t
contains
procedure(get_scalar_field), deferred :: get_scalar_field
end type scalar_field_generator_t

type, public, abstract :: vector_field_generator_t
contains
procedure(get_vector_field), deferred :: get_vector_field
end type vector_field_generator_t

!!!!!!!!!!!!Specific fields
type, extends(scalar_field_generator_t) :: xyz_scalar_field_generator_t
contains
procedure :: get_scalar_field => generate_xyz_scalar_field
end type xyz_scalar_field_generator_t

type, extends(vector_field_generator_t) :: solid_rotation_field_generator_t
contains
procedure :: get_vector_field => generate_solid_rotation_vector_field
end type solid_rotation_field_generator_t

type, extends(vector_field_generator_t) :: xyz_grad_generator_t
contains
procedure :: get_vector_field => generate_xyz_grad_field
end type xyz_grad_generator_t

!!!field generator instances
type(xyz_scalar_field_generator_t)     :: xyz_scalar_field_generator
type(solid_rotation_field_generator_t) :: solid_rotation_field_generator
type(xyz_grad_generator_t)             :: xyz_grad_generator

abstract interface
    subroutine get_scalar_field(this,f,npts,nlev,x,y,z)
        import scalar_field_generator_t
        class(scalar_field_generator_t),  intent(in) :: this
        integer(kind=4), intent(in)                  :: npts, nlev
        real(kind=8),    intent(in)                  :: x(npts), y(npts), z(npts)
        real(kind=8),    intent(out)                 :: f(npts,nlev)
    end subroutine get_scalar_field
    subroutine get_vector_field(this,vx,vy,vz,npts,nlev,x,y,z)
        import vector_field_generator_t
        class(vector_field_generator_t),    intent(in)  :: this
        integer(kind=4),                    intent(in)  :: npts, nlev
        real(kind=8),                       intent(in)  :: x(npts), y(npts), z(npts)
        real(kind=8), dimension(npts,nlev), intent(out) :: vx, vy, vz
    end subroutine get_vector_field
end interface

contains

subroutine set_scalar_test_field(f, generator, mesh, halo_width, fill_value)
    type(grid_field_t),              intent(inout) :: f
    class(scalar_field_generator_t), intent(in)    :: generator
    type(mesh_t),                    intent(in)    :: mesh
    integer(kind=4),                 intent(in)    :: halo_width
    real(kind=8),          optional, intent(in)    :: fill_value 

    !locals
    integer(kind=4) t, isv, iev, jsv, jev, npoints, klev
    real(kind=8), allocatable :: p(:,:,:)
    real(kind=8), dimension(:,:), allocatable :: px, py, pz 

    do t = mesh%ts, mesh%te
        isv = mesh%tile(t)%is-halo_width
        iev = mesh%tile(t)%ie+halo_width
        jsv = mesh%tile(t)%js-halo_width
        jev = mesh%tile(t)%je+halo_width
        npoints = (iev-isv+1)*(jev-jsv+1)
        klev = mesh%tile(t)%ke-mesh%tile(t)%ks+1
        if(present(fill_value)) f%tile(t)%p(:,:,:) = -fill_value
        !suppress temporary array warning (made tmp arrays manually)
        allocate(p(isv:iev,jsv:jev,1:klev))
        allocate(px(isv:iev,jsv:jev),py(isv:iev,jsv:jev),pz(isv:iev,jsv:jev))
        px(isv:iev,jsv:jev) = mesh%tile(t)%rx(isv:iev,jsv:jev)
        py(isv:iev,jsv:jev) = mesh%tile(t)%ry(isv:iev,jsv:jev)
        pz(isv:iev,jsv:jev) = mesh%tile(t)%rz(isv:iev,jsv:jev)
        call generator%get_scalar_field(p, npoints,klev,px,py,pz)
        f%tile(t)%p(isv:iev,jsv:jev,1:klev) = p(isv:iev,jsv:jev,1:klev)
        deallocate(p,px,py,pz)
    end do
end subroutine set_scalar_test_field

subroutine set_vector_test_field(u, v, generator, mesh_u, mesh_v, halo_width, components_type, fill_value)

    use parcomm_mod, only : parcomm_global

    type(grid_field_t),              intent(inout) :: u, v
    class(vector_field_generator_t), intent(in)    :: generator    
    type(mesh_t),                    intent(in)    :: mesh_u, mesh_v
    integer(kind=4),                 intent(in)    :: halo_width
    character(len=*),                intent(in)    :: components_type
    real(kind=8),          optional, intent(in)    :: fill_value
    
    !locals
    integer(kind=4) :: t

     do t = mesh_u%ts, mesh_u%te
        if(components_type=="contravariant") then
            call set_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%b1, &
                                                   halo_width,fill_value)
            call set_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%b2, &
                                                   halo_width,fill_value)
        else if (components_type=="covariant") then
            call set_vector_test_field_1tile_1comp(u%tile(t),generator,mesh_u%tile(t),mesh_u%tile(t)%a1, &
                                                   halo_width,fill_value)
            call set_vector_test_field_1tile_1comp(v%tile(t),generator,mesh_v%tile(t),mesh_v%tile(t)%a2, &
                                                   halo_width,fill_value)
        else
            call parcomm_global%abort("test_functions_mod, set_vector_test_field: unknown components type "// &
                                      components_type//" use covariant or contravariant")
        end if

     end do

end subroutine set_vector_test_field

subroutine set_vector_test_field_1tile_1comp(u,generator,mesh,bvec,halo_width,fill_value)
    use grid_field_mod, only : tile_field_t
    use mesh_mod,       only : tile_mesh_t

    type(tile_field_t),              intent(inout) :: u
    class(vector_field_generator_t), intent(in)    :: generator    
    type(tile_mesh_t),               intent(in)    :: mesh
    real(kind=8),                    intent(in)    :: bvec(1:3,mesh%is-mesh%halo_width:mesh%ie+mesh%halo_width, &
                                                               mesh%js-mesh%halo_width:mesh%je+mesh%halo_width)
    integer(kind=4),                 intent(in)    :: halo_width
    real(kind=8),          optional, intent(in)    :: fill_value
    
    !locals
    integer(kind=4) :: isv, iev, jsv, jev, i, j, k, k1, klev, npoints
    real(kind=8), dimension(:,:,:), allocatable ::  vx, vy, vz
    real(kind=8), dimension(:,:), allocatable   :: px, py, pz

    if(present(fill_value)) u%p(:,:,:) = fill_value
    
    isv = mesh%is-halo_width
    iev = mesh%ie+halo_width
    jsv = mesh%js-halo_width
    jev = mesh%je+halo_width
    klev = mesh%ke-mesh%ks+1

    allocate(vx(isv:iev,jsv:jev,1:klev),vy(isv:iev,jsv:jev,1:klev),vz(isv:iev,jsv:jev,1:klev))
    allocate(px(isv:iev,jsv:jev),py(isv:iev,jsv:jev),pz(isv:iev,jsv:jev))
    px(isv:iev,jsv:jev) = mesh%rx(isv:iev,jsv:jev)
    py(isv:iev,jsv:jev) = mesh%ry(isv:iev,jsv:jev)
    pz(isv:iev,jsv:jev) = mesh%rz(isv:iev,jsv:jev)
    npoints = (iev-isv+1)*(jev-jsv+1)
    call generator%get_vector_field(vx,vy,vz,npoints,klev,px,py,pz)

    do k=mesh%ks, mesh%ke
        k1 = k-mesh%ks+1
        do j = jsv, jev
            do i = isv, iev
            u%p(i,j,k) = sum([vx(i,j,k1),vy(i,j,k1),vz(i,j,k1)]*bvec(1:3,i,j))
            end do
        end do
    end do
    deallocate(vx,vy,vz)
    deallocate(px,py,pz)
end subroutine set_vector_test_field_1tile_1comp

subroutine generate_xyz_scalar_field(this,f,npts,nlev,x,y,z)
    import scalar_field_generator_t
    class(xyz_scalar_field_generator_t),  intent(in) :: this
    integer(kind=4), intent(in)                  :: npts, nlev
    real(kind=8),    intent(in)                  :: x(npts), y(npts), z(npts)
    real(kind=8),    intent(out)                 :: f(npts,nlev)

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: M(3,3) = reshape([1.0_8, 0.0_8, 0.0_8, &
                                        0.0_8, 1.0_8, 0.0_8,  &
                                        0.0_8, 0.0_8, 1.0_8],[3,3])

    do k = 1, nlev
        k1 = mod(k-1,3)+1 !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            f(i,k) = M(1,k1)*x(i)+M(2,k1)*y(i)+M(3,k1)*z(i) !k1=1:x, k1=2:y, k1=3:z
        end do
    end do
end subroutine generate_xyz_scalar_field

subroutine generate_solid_rotation_vector_field(this,vx,vy,vz,npts,nlev,x,y,z)
    class(solid_rotation_field_generator_t),  intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: axis(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & !rotation axes for different levs
                                            0.0_8, 1.0_8, 0.0_8,  &
                                            0.0_8, 0.0_8, 1.0_8],[3,3])

    do k = 1, nlev
        k1 = mod(k-1,3)+1  !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            vx(i,k) = axis(2,k1)*z(i)-axis(3,k1)*y(i)
            vy(i,k) =-axis(1,k1)*z(i)+axis(3,k1)*x(i)
            vz(i,k) = axis(1,k1)*y(i)-axis(2,k1)*x(i)
        end do
    end do
end subroutine generate_solid_rotation_vector_field

subroutine generate_xyz_grad_field(this,vx,vy,vz,npts,nlev,x,y,z)
    class(xyz_grad_generator_t),              intent(in)  :: this
    integer(kind=4),                          intent(in)  :: npts, nlev
    real(kind=8),       dimension(npts),      intent(in)  :: x, y, z
    real(kind=8),       dimension(npts,nlev), intent(out) :: vx, vy, vz

    integer(kind=4) :: i, k, k1
    real(kind=8)    :: e(3,3) = reshape([1.0_8, 0.0_8, 0.0_8,  & ! level1: f=x, 
                                         0.0_8, 1.0_8, 0.0_8,  & ! level2: f=y
                                         0.0_8, 0.0_8, 1.0_8],[3,3]) !level 3: f=z
    !(nabla)_h f = (nabla)_3d f- n *(n, nabla_3d f)
    !(nabla)_3d f = e(k)
    real(kind=8)    :: n(3) !normal to the surface
    real(kind=8)    :: nne(3) !n *(n, nabla_3d f)
    real(kind=8)    :: nabla_h(3)

    do k = 1, nlev
        k1 = mod(k-1,3)+1  !1,2,3,1,2,3,1,2 etc
        do i=1, npts
            n(1:3) = [x(i),y(i),z(i)] / sqrt(x(i)**2+y(i)**2+z(i)**2)
            nne(1:3) = n(1:3)*sum(n(1:3)*e(1:3,k))
            nabla_h(1:3) = e(1:3,k)-nne(1:3)
            vx(i,k) = nabla_h(1)
            vy(i,k) = nabla_h(2)
            vz(i,k) = nabla_h(3)
        end do
    end do
end subroutine generate_xyz_grad_field

end module test_fields_mod
