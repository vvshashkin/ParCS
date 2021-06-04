module metric_mod

implicit none

!Abstract type to aquire horizontal grid characteristics from panel coordinates
type, abstract :: metric_t
    real(kind=8) a      ! projection scale factor (radius of the sphere in the case of equiangular cubed sphere)
    real(kind=8) x0, y0 ! lower bound of panel coordinates
    real(kind=8) x1, y1 ! upper bound of panel coordinates

contains

    procedure(vector_cart), deferred :: point_r !Cartesian 3d coordinates of points
    !3d cartesian coordinates of grid basis vectors
    procedure(vector_cart), deferred :: a1 !covariant dr/dx
    procedure(vector_cart), deferred :: a2 !covariant dr/dy
    procedure(vector_cart), deferred :: b1 !contravariant x direction
    procedure(vector_cart), deferred :: b2 !contravariant y direction
    !Metric tensor etc
    procedure(symtensor2),  deferred :: Q !metric tensor (a1*a1 a1*a2; a1*a2 a2*a2)
    procedure(symtensor2),  deferred :: QI !inversed metric tensor
    procedure(tensor0),     deferred :: G !sqrt of metric tensor det

end type metric_t

abstract interface
    function vector_cart(this,panel_ind, x, y) result(r)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: x, y
        real(kind=8)                :: r(3)
    end function vector_cart

    function symtensor2(this,panel_ind,x,y) result(Q)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: x, y
        real(kind=8)                :: Q(3)
    end

    function tensor0(this,panel_ind,x,y) result(G)
        import metric_t
        class(metric_t), intent(in) :: this
        integer(kind=4), intent(in) :: panel_ind
        real(kind=8),    intent(in) :: x, y
        real(kind=8)                :: G
    end
end interface

end module metric_mod
