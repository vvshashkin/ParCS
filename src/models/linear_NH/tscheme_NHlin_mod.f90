module tscheme_NHlin_mod

implicit none

contains

subroutine init_tscheme_NHlin(time_scheme, operator, model_params, stvec, &
                              master_id, myid, np, namelist_str)
    use timescheme_abstract_mod, only : timescheme_abstract_t
    use operator_NHlin_mod,      only : operator_NHlin_t
    use parameters_NHlin_mod,    only : parameters_NHlin_t
    use stvec_NHlin_mod,         only : stvec_NHlin_t
    !time schemes
    use rk4_mod,                  only : init_rk4
    use exp_krylov_mod,           only : init_exp_krylov
    use ars232_mod,               only : init_ars232

    class(timescheme_abstract_t), allocatable, &
                                  intent (inout) :: time_scheme
    class(operator_NHlin_t),      intent(in)     :: operator
    class(parameters_NHlin_t),    intent(in)     :: model_params
    class(stvec_NHlin_t),         intent(in)     :: stvec
    integer(kind=4),              intent(in)     :: master_id, myid, np
    character(:), allocatable,    intent(in)     :: namelist_str

    character(256)  :: time_scheme_name = "RK4"
    integer(kind=4) :: Mmax=30, iom=2 ! for Krylov exp scheme

    namelist /tscheme/ time_scheme_name, Mmax, iom

    if(allocated(namelist_str)) then
        read(namelist_str,tscheme)
    end if

    if(trim(time_scheme_name) == "RK4") then
        print *, "Using RK4 time stepping scheme"
        time_scheme = init_rk4(operator, stvec)
    else if(trim(time_scheme_name) == "ARS232") then
        print *, "Using ARS232 time stepping scheme"
        call init_ars232(time_scheme, operator, stvec)
    else if(trim(time_scheme_name) == "EXP") then
        print '(A,I5,A,I4)', "Using Krylov exponential time stepping scheme, Mmax=", Mmax, ", iom=", iom
        call init_exp_krylov(time_scheme, operator, stvec, 5)
    else
        call avost("init_tscheme_NHlin: time scheme " // trim(time_scheme_name) // &
                                             " is unknown or not implemented")
    end if

end subroutine init_tscheme_NHlin

end module tscheme_NHlin_mod
