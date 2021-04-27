module swlin_operator_factory_mod

implicit none

character(64) :: div_op_name = "divergence2", grad_op_name = "gradient2"
namelist /oper_ini/ div_op_name, grad_op_name

contains

subroutine create_swlin_operator(swlin_operator, model_params, master_id, myid, np, namelist_str)

    use partition_mod,        only : partition_t
    use exchange_factory_mod, only : create_Agrid_halo_exchange
    use ecs_halo_factory_mod, only : init_ecs_halo
    use hor_difops_basic_mod, only : cl_gradient_contra_c2, cl_divergence_cgr2, &
                                     cl_gradient_0, cl_divergence_0
    use parameters_swlin_mod, only : parameters_swlin_t
    use operator_swlin_mod,   only : operator_swlin_t

    type(operator_swlin_t)                 :: swlin_operator
    type(parameters_swlin_t),  intent(in)  :: model_params
    integer(kind=4),           intent(in)  :: master_id, myid, np
    character(:), allocatable, intent(in)  :: namelist_str

    integer(kind=4) ind

    swlin_operator%exch_halo = create_Agrid_halo_exchange(model_params%partition,    &
                                                model_params%halo_width,   &
                                                'full', myid, np)
    allocate(swlin_operator%halo(model_params%ts:model_params%te))
    do ind = model_params%ts, model_params%te
        swlin_operator%halo(ind) = init_ecs_halo(model_params%mesh(ind)%is, model_params%mesh(ind)%ie, &
                                       model_params%mesh(ind)%js, model_params%mesh(ind)%je, &
                                       model_params%mesh(ind)%nx, swlin_operator%op_halo_width,&
                                       model_params%mesh(ind)%hx)
    end do

    if(allocated(namelist_str)) then
        read(namelist_str, oper_ini)
    end if

    if(trim(grad_op_name) == "gradient2") then
        swlin_operator%grad_contra => cl_gradient_contra_c2
    else if(trim(grad_op_name) == "gradient0") then
        swlin_operator%grad_contra => cl_gradient_0
    else
        call avost("SWLIN operator init: unknown gradient operator -- " // &
                                                          trim(grad_op_name))
    end if

    if(trim(div_op_name) == "divergence2") then
        swlin_operator%div         => cl_divergence_cgr2
    else if(trim(div_op_name) == "divergence0") then
        swlin_operator%div         => cl_divergence_0
    else
        call avost("SWLIN operator init: unknown gradient operator -- " // &
                                                          trim(div_op_name))
    end if

    if(myid == master_id) then
        print *, "---Operator initialization"
        print *, "gradient operator: ", grad_op_name
        print *, "divergence operator: ", div_op_name
        print *, "--------------------------"
    end if

end subroutine create_swlin_operator

end module swlin_operator_factory_mod
