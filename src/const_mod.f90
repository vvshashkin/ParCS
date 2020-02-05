!Module for fundamental constants
module const_mod
implicit none
real(kind=8), parameter :: pi = acos(-1._8)
real(kind=8), parameter :: grav = 9.80616_8   !gravity acceleration m/s^2
real(kind=8), parameter :: radz = 0.6371229d7 !Earth radius
real(kind=8), parameter :: Day24h_sec = 24._8*3600._8 !day length in seconds
real(kind=8)               gs_axis(3)         !geographical south pole to north pole axis (of unit length)

contains
subroutine const_mod_init()
   gs_axis = [0._8, 0._8, 1._8]
end subroutine const_mod_init

end module const_mod
