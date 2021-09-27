module test_latlon_output_mod

use domain_mod,             only : domain_t
use domain_factory_mod,     only : create_domain
use grid_field_mod,         only : grid_field_t
use grid_field_factory_mod, only : create_grid_field
use outputer_abstract_mod,  only : outputer_t
use outputer_factory_mod,   only : create_latlon_outputer
use parcomm_mod,            only : parcomm_global

implicit none

contains

subroutine test_latlon_output(staggering, scalar_grid)

    use test_fields_mod,        only : set_scalar_test_field, &
                                       xyz_f=>xyz_scalar_field_generator

    character(len=*), intent(in) :: staggering, scalar_grid

    class(outputer_t), allocatable  :: outputer
    type(domain_t)                  :: domain
    type(grid_field_t)              :: f, u, v

    character(*), parameter    :: file_name = "h.dat"
    integer(kind=4), parameter :: nh = 100, nz = 3

    call create_domain(domain, "cube", staggering, nh, nz)

    if(scalar_grid == "A") then
        call create_grid_field(f, 0, 0, domain%mesh_o)
        call set_scalar_test_field(f,xyz_f,domain%mesh_o,0)
    else if(scalar_grid == "Ah") then
        call create_grid_field(f, 0, 0, domain%mesh_xy)
        call set_scalar_test_field(f,xyz_f,domain%mesh_xy,0)
    else
        call parcomm_global%abort("test_latlon_outputter_mod, "//&
                                  "unsupported type of scalar point:"//&
                                  scalar_grid)
    end if

    call create_grid_field(u, 0, 0, domain%mesh_u)
    call create_grid_field(v, 0, 0, domain%mesh_v)

    call create_latlon_outputer(outputer, 201, 400, scalar_grid, domain)
    call outputer%write(f,domain,"h.dat",1)

    if(domain%parcomm%myid == 0) then
        print *, "Test latlon outputer, staggering: "//staggering //&
                " scalar  grid: "// scalar_grid//" passed"
    end if

end subroutine test_latlon_output

end module test_latlon_output_mod
