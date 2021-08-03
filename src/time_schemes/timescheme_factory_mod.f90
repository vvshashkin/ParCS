module timescheme_factory_mod

use stvec_mod,         only : stvec_t
use timescheme_mod,    only : timescheme_t
use parcomm_mod,       only : parcomm_global
use explicit_Eul1_mod, only : init_explicit_Eul1_ts
use rk4_mod,           only : init_rk4

implicit none

contains

subroutine create_timescheme(tscheme, v, tscheme_name)
    class(timescheme_t), allocatable, intent(out) :: tscheme
    class(stvec_t),   intent(in) :: v !example of model state-vector
    character(len=*), intent(in) :: tscheme_name

    if(tscheme_name == "explicit_Eul1") then
        call init_explicit_Eul1_ts(tscheme,v)
    elseif(tscheme_name == "rk4") then
        call init_rk4(tscheme,v)
    else
        call parcomm_global%abort("unknown timescheme_name in create_timescheme: "// &
                                  tscheme_name)
    end if

end subroutine create_timescheme

end module timescheme_factory_mod
