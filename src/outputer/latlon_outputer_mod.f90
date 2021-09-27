module latlon_outputer_mod

use outputer_abstract_mod, only : outputer_t
use grid_field_mod,        only : grid_field_t
use exchange_abstract_mod, only : exchange_t
use domain_mod,            only : domain_t
use tiles_mod,             only : tiles_t
use abstract_regrid_mod,   only : regrid_t

implicit none

type, public, extends(outputer_t) :: latlon_outputer_t
    class(exchange_t), allocatable :: gather_exch
    type(grid_field_t)             :: exchange_buf
    integer(kind=4)                :: master_id
    type(tiles_t)                  :: tiles
    class(regrid_t),   allocatable :: regrid
    type(domain_t)                 :: regrid_domain
    type(grid_field_t)             :: regrid_work
    integer(kind=4)                :: Nlon, Nlat
    real(kind=8),      allocatable :: buffer(:,:,:)

contains
    procedure, public :: write => latlon_write
end type latlon_outputer_t

contains

subroutine latlon_write(this, f, domain, file_name, rec_num)

    class(latlon_outputer_t),         intent(inout) :: this
    type(grid_field_t),               intent(inout) :: f
    character(*),                     intent(in)    :: file_name
    type(domain_t),                   intent(in)    :: domain
    integer(kind=4),  optional,       intent(in)    :: rec_num

    integer(kind=4) :: i, j, k, t, ks, ke, panel_ind
    integer(kind=4) :: ts, te, js, je, is, ie

    if (domain%parcomm%myid == this%master_id) then
        do t = domain%partition%ts, domain%partition%te
            do k = this%tiles%ks(t), this%tiles%ke(t)
                do j = this%tiles%js(t), this%tiles%je(t)
                    do i = this%tiles%is(t), this%tiles%ie(t)
                        this%exchange_buf%tile(t)%p(i,j,k) = f%tile(t)%p(i,j,k)
                    end do
                end do
            end do
        end do
        call this%gather_exch%do(this%exchange_buf, domain%parcomm)
    else
        call this%gather_exch%do(f, domain%parcomm)
    end if

    if(domain%parcomm%myid == this%master_id) then

        open(newunit=this%out_stream, file = trim(file_name), &
             access="direct", recl = this%Nlat*this%Nlon)

        ts = this%exchange_buf%ts
        te = this%exchange_buf%te
        ks = this%exchange_buf%tile(ts)%ks
        ke = this%exchange_buf%tile(ts)%ke
        !map exchanged grid function to work buffer before regrid
        do k = ks,ke

            do t=ts,te
                is = this%exchange_buf%tile(t)%is
                ie = this%exchange_buf%tile(t)%ie
                js = this%exchange_buf%tile(t)%js
                je = this%exchange_buf%tile(t)%je
                panel_ind = domain%partition%panel_map(t)
                do j=js, je; do i=is, ie
                    this%regrid_work%tile(panel_ind)%p(i,j,1) = this%exchange_buf%tile(t)%p(i,j,k)
                end do; end do
            end do

            call this%regrid%do_regrid(this%buffer,this%regrid_work,this%regrid_domain)
            write(this%out_stream,rec=(rec_num-1)*(ks-ke)+(k-ks+1)) real(this%buffer,4)
            !print *, k, size(this%buffer), this%Nlat*this%Nlon
        end do
        close(this%out_stream)
    end if

end subroutine latlon_write

end module latlon_outputer_mod
