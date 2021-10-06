module quadrature_factory_mod

use abstract_quadrature_mod, only : quadrature_t
use default_quadrature_mod,  only : default_tile_quadrature_t
use parcomm_mod,             only : parcomm_global
use domain_mod,              only : domain_t

implicit none

contains

subroutine create_quadrature(quadrature, quadrature_name, domain)

    class(quadrature_t), allocatable, intent(out) :: quadrature
    character(len=*),                 intent(in)  :: quadrature_name
    type(domain_t),                   intent(in)  :: domain

    integer(kind=4) :: ts, te

    allocate(quadrature)
    if(quadrature_name == "default_quadrature") then
        ts = domain%partition%ts
        te = domain%partition%te
        allocate(default_tile_quadrature_t :: quadrature%tile(ts:te))
    else if(quadrature_name == "SBP_Ah21_quadrature" .or. &
            quadrature_name == "SBP_Ah42_quadrature") then
        call create_ah_sbp_quadrature(quadrature, quadrature_name, domain)
    ! else if(quadrature_name == "SBP_C21_quadrature" .or. &
    !         quadrature_name == "SBP_C42_quadrature") then
    !     call create_C_sbp_quadrature(quadrature, quadrature_name, domain)
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_quadrature, quadrature_factory_mod")
    end if

end subroutine create_quadrature

subroutine create_ah_sbp_quadrature(quadrature, quadrature_name, domain)

    use sbp_quadrature_mod, only : sbp_tile_quadrature_t
    use sbp_operators_collection_mod, only : Q21_A, Q42_A

    class(quadrature_t), allocatable, intent(inout) :: quadrature
    character(len=*),                 intent(in)    :: quadrature_name
    type(domain_t),                   intent(in)    :: domain

    type(sbp_tile_quadrature_t), allocatable :: tile_q(:)
    integer(kind=4) :: t, is, ie, js, je, i, j, n
    real(kind=8) A_edge(8)
    integer(kind=4) :: nedge

    if(quadrature_name == "SBP_Ah21_quadrature") then
        nedge = size(Q21_A)
        A_edge(1:nedge) =  Q21_A(1:nedge)
    else if(quadrature_name == "SBP_Ah42_quadrature") then
        nedge = size(Q42_A)
        A_edge(1:nedge) =  Q42_A(1:nedge)
    else
        call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
                                  " in create_sbp_quadrature, quadrature_factory_mod")
    end if

    allocate(tile_q(domain%partition%ts:domain%partition%te))

    do t = domain%partition%ts, domain%partition%te
       is = domain%mesh_xy%tile(t)%is; ie = domain%mesh_xy%tile(t)%ie
       js = domain%mesh_xy%tile(t)%js; je = domain%mesh_xy%tile(t)%je

       allocate(tile_q(t)%Ax(is:ie))
       allocate(tile_q(t)%Ay(js:je))

       call create_mass_coefficients(tile_q(t)%Ax, is, ie, domain%mesh_xy%tile(t)%nx, &
                                                                 A_edge(1:nedge), nedge)
       call create_mass_coefficients(tile_q(t)%Ay, js, je, domain%mesh_xy%tile(t)%ny, &
                                                                 A_edge(1:nedge), nedge)
   end do

   call move_alloc(tile_q, quadrature%tile)
end

! subroutine create_C_sbp_quadrature(quadrature, quadrature_name, domain)
!
!     use sbp_quadrature_mod, only : sbp_tile_quadrature_t
!     use sbp_operators_collection_mod, only : D42_A_interfaces, D42_A_centers, &
!                                              D21_A_interfaces, D21_A_centers
!
!     class(quadrature_t), allocatable, intent(inout) :: quadrature
!     character(len=*),                 intent(in)    :: quadrature_name
!     type(domain_t),                   intent(in)    :: domain
!
!     type(sbp_tile_quadrature_t), allocatable :: tile_q(:)
!     integer(kind=4) :: t, is, ie, js, je, i, j, n
!     real(kind=8) A_edge(8)
!     integer(kind=4) :: nedge
!
!     if(quadrature_name == "SBP_Ah21_quadrature") then
!         A_edge(1:2) =  Q21_A(1:2)
!         nedge = 2
!     else if(quadrature_name == "SBP_Ah42_quadrature") then
!         A_edge(1:4) =  Q42_A(1:4)
!         nedge = 4
!     else
!         call parcomm_global%abort("unknown quadrature:"// quadrature_name //&
!                                   " in create_sbp_quadrature, quadrature_factory_mod")
!     end if
!
!     allocate(tile_q(domain%partition%ts:domain%partition%te))
!
!     do t = domain%partition%ts, domain%partition%te
!        is = domain%mesh_xy%tile(t)%is; ie = domain%mesh_xy%tile(t)%ie
!        js = domain%mesh_xy%tile(t)%js; je = domain%mesh_xy%tile(t)%je
!
!        allocate(tile_q(t)%Ax(is:ie))
!        allocate(tile_q(t)%Ay(js:je))
!
!        call create_mass_coefficients(tile_q(t)%Ax, is, ie, domain%mesh_xy%tile(t)%nx, &
!                                                                  A_edge(1:nedge), nedge)
!        call create_mass_coefficients(tile_q(t)%Ay, js, je, domain%mesh_xy%tile(t)%ny, &
!                                                                  A_edge(1:nedge), nedge)
!    end do
!
!    call move_alloc(tile_q, quadrature%tile)
! end


subroutine create_mass_coefficients(A, is, ie, n, A_edge, n_edge)

    integer(kind=4), intent(in)    :: is, ie, n, n_edge
    real(kind=8),    intent(in)    :: A_edge(n_edge)

    real(kind=8),    intent(inout) :: A(is:ie)

    integer(kind=8) :: i

    A(is:ie) = 1.0_8

    do i = max(1, is), min(n_edge,ie)
        A(i) = A_edge(i)
    end do

    do i = max(n-n_edge+1, is), min(n,ie)
        A(i) = A_edge(n-i+1)
    end do

end subroutine create_mass_coefficients

end module quadrature_factory_mod
