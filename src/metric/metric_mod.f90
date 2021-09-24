module metric_mod

implicit none

!Abstract type to aquire horizontal grid characteristics from panel coordinates
type, abstract :: metric_t
    real(kind=8) scale      ! projection scale factor (radius of the sphere in the case of equiangular cubed sphere)
    real(kind=8) alpha0, beta0 ! lower bound of panel coordinates
    real(kind=8) alpha1, beta1 ! upper bound of panel coordinates

    real(kind=8) :: omega ! rotation speed of the coord system
    real(kind=8), allocatable :: rotation_axis(:)!direction of the angular velocity vector

    real(kind=8), allocatable :: rotation_matrix(:,:)

contains

    procedure(vector_cart),    deferred :: point_r !Cartesian 3d coordinates of points
    procedure(cart_to_native_transform), deferred :: transform_cartesian_to_native
    !3d cartesian coordinates of grid basis vectors
    procedure(vector_cart),    deferred :: a1 !covariant dr/dx
    procedure(vector_cart),    deferred :: a2 !covariant dr/dy
    procedure(vector_cart),    deferred :: b1 !contravariant x direction
    procedure(vector_cart),    deferred :: b2 !contravariant y direction
    !Metric tensor etc
    procedure(symtensor2),     deferred :: Q !metric tensor (a1*a1 a1*a2; a1*a2 a2*a2)
    procedure(symtensor2),     deferred :: QI !inversed metric tensor
    procedure(tensor0),        deferred :: G !sqrt of metric tensor det

end type metric_t

abstract interface
    function vector_cart(this,panel_ind, alpha, beta) result(r)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: r(3)
    end function vector_cart

    subroutine cart_to_native_transform(this,panel_ind, alpha, beta, r)
        import metric_t
        class(metric_t), intent(in)  :: this
        integer(kind=4), intent(out) :: panel_ind
        real(kind=8),    intent(out) :: alpha, beta
        real(kind=8),    intent(in)  :: r(3)
    end subroutine cart_to_native_transform

    function symtensor2(this,panel_ind,alpha,beta) result(Q)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: Q(3)
    end

    function tensor0(this,panel_ind,alpha,beta) result(G)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: alpha, beta
        real(kind=8)                :: G
    end
end interface

end module metric_mod
