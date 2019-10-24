module partition_mod
use tile_mod, only : tile_t
implicit none

type, public :: partition_t
    type(tile_t),    allocatable :: tile(:)     !array of partition tiles
    integer(kind=4), allocatable :: proc_map(:) !determine belonging of the tile to the specific processor
    integer(kind=4)              :: npoints     !number of grid points in x/y direction for the one panel
    integer(kind=4)              :: num_tiles   !number of tiles in the partition
contains
    procedure, public :: init
    procedure, public :: write_to_txt
end type partition_t

contains

subroutine init(this, Nh, Nz, num_tiles, strategy)
    class(partition_t), intent(inout) :: this
    integer(kind=4),    intent(in)    :: Nh        ! num of points in x and y direction at each panel
    integer(kind=4),    intent(in)    :: Nz        ! num of points in z direction at each panel
    integer(kind=4),    intent(in)    :: num_tiles ! num of tiles at each panel
    character(*),       intent(in)    :: strategy

!We assume that number of tiles at each panel are the same

    allocate(this%tile(6*num_tiles))
    allocate(this%proc_map(6*num_tiles))
    this%num_tiles = num_tiles
    this%npoints = Nh

    select case (strategy)

    case ('default')
        call default_strategy(this, Nh, Nz)

    case default
        print*, 'Wrong strategy in partition_mod.f90'
        stop

    end select
end subroutine init

subroutine default_strategy(partition, Nh, Nz)
    type(partition_t), intent(inout) :: partition
    integer(kind=4),    intent(in)   :: Nh
    integer(kind=4),    intent(in)   :: Nz

    integer(kind=4) :: npx_npy(2)
    integer(kind=4), allocatable :: wx(:), wy(:)
    integer(kind=4) :: ind, i, j, panel_ind
    integer(kind=4) :: is, ie, js, je, ks, ke

    npx_npy = get_factors(partition%num_tiles)

    allocate(wx(npx_npy(1)), wy(npx_npy(2)))

    call partition_1d(npx_npy(1), Nh, wx)
    call partition_1d(npx_npy(2), Nh, wy)

    ind = 0
    do panel_ind = 1, 6
        do j = 1, npx_npy(2)
            do i = 1, npx_npy(1)

                ind = ind+1

                if (j==1) then
                    js = 1
                    je = wy(j)
                else
                    js = partition%tile(ind-npx_npy(1))%js + wy(j-1)
                    je = partition%tile(ind-npx_npy(1))%je + wy(j  )
                end if

                if (i==1) then
                    is = 1
                    ie = wx(i)
                else
                    is = partition%tile(ind-1)%is + wx(i-1)
                    ie = partition%tile(ind-1)%ie + wx(i  )
                endif

                ks = 1
                ke = Nz

                call partition%tile(ind)%init(is, ie, js, je, ks, ke, panel_ind)

                ! call partition%tile(ind)%print()

                partition%proc_map(ind) = ind-1

            end do
        end do
    end do

    call partition%write_to_txt('partition.txt')

end subroutine default_strategy

subroutine write_to_txt(this, filename)
    class(partition_t), intent(in) :: this
    character(*),       intent(in) :: filename

    integer(kind=4) :: ind

    open(120, file=filename)

    do ind = 1, size(this%tile,1)
        write(120,*) this%tile(ind)%panel_number, this%tile(ind)%is, this%tile(ind)%ie &
                                           , this%tile(ind)%js, this%tile(ind)%je, this%proc_map(ind)
    end do

    close(120)

end subroutine

subroutine partition_1d(Np, N, work)
    integer(kind=4), intent(in)  :: Np, N
    integer(kind=4), intent(out) :: work(Np)

    real(kind=8) :: mean_work, deficite

    mean_work = N/Np
    work(1:Np)=mean_work
    deficite = mod(N,Np)
    work(1:deficite) = work(1:deficite)+1

end subroutine partition_1d

function get_factors(n) result(factors)
    integer(kind=4), intent(in)  :: n
    integer(kind=4) :: factors(2)

    integer(kind=4) :: f1

    f1 = int(sqrt(real(n,8)))

    do
        if (mod(n,f1)==0) exit
        f1 = f1-1
    end do

    factors = (/n/f1, f1/)

end function get_factors

end module partition_mod
