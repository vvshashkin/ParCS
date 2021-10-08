module key_value_mod

implicit none

type string_t
    character(:), allocatable :: str
end type string_t

type key_value_r8_t
    type(string_t), allocatable :: keys(:)
    real(kind=8),   allocatable :: values(:)
end type key_value_r8_t

end module key_value_mod
