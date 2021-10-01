module sbp_factory_mod

use sbp_operator_mod,             only : sbp_operator_t
use parcomm_mod,                  only : parcomm_global

use sbp_operators_collection_mod, only : Q21, lastnonzeroQ21, Da2_in, Da2_inshift, &
                                         Q42, lastnonzeroQ42, Da4_in, Da4_inshift, &
                                         Q43, lastnonzeroQ43,                      &
                                         W42_staggered_i2c_noopt, W42_staggered_i2c_opt1, &
                                         W42_staggered_i2c_opt2, W42_staggered_i2c_last_nonzero, &
                                         W42_staggered_c2i_noopt, W42_staggered_c2i_opt1, &
                                         W42_staggered_c2i_opt2, W42_staggered_c2i_last_nonzero, &
                                         W42_staggered_in, W42_staggered_c2i_in_shift, &
                                         W42_staggered_c2i_in_shift,  W42_staggered_i2c_in_shift

implicit none

contains

function create_sbp_operator(sbp_operator_name) result(sbp_op)
    character(len=*), intent(in) :: sbp_operator_name
    type(sbp_operator_t)         :: sbp_op

    !the sign of right corner part of matrix
    !should be -1 for derivatives, 1 for interpolations
    real(kind=8)                 :: right_side_sign
    !size of edge matrix block
    integer(kind=4)              :: n_edge, nw_edge
    !

    !Initialization of Left-side and inner matrix coefficients
    !Right side matrix block is initialized below for all options
    !using symmetry
    if(sbp_operator_name == "d21") then
        sbp_op%W_edge_l    = Q21
        sbp_op%edge_last_l = lastnonzeroQ21
        sbp_op%W_in        = Da2_in
        sbp_op%in_shift    = Da2_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "d42") then
        sbp_op%W_edge_l    = Q42
        sbp_op%edge_last_l = lastnonzeroQ42
        sbp_op%W_in        = Da4_in
        sbp_op%in_shift    = Da4_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "d43") then
        sbp_op%W_edge_l    = Q43
        sbp_op%edge_last_l = lastnonzeroQ43
        sbp_op%W_in        = Da4_in
        sbp_op%in_shift    = Da4_inshift
        sbp_op%dnx = 0
        right_side_sign = -1.0_8
    else if(sbp_operator_name == "W42_stagered_interp_c2i") then
        sbp_op%W_edge_l    = W42_staggered_c2i_opt1
        sbp_op%edge_last_l = W42_staggered_c2i_last_nonzero
        sbp_op%W_in        = W42_staggered_in
        sbp_op%in_shift    = W42_staggered_c2i_in_shift
        sbp_op%dnx = 1
        right_side_sign = 1.0_8
    else if(sbp_operator_name == "W42_stagered_interp_i2c") then
        sbp_op%W_edge_l    = W42_staggered_i2c_opt1
        sbp_op%edge_last_l = W42_staggered_i2c_last_nonzero
        sbp_op%W_in        = W42_staggered_in
        sbp_op%in_shift    = W42_staggered_i2c_in_shift
        sbp_op%dnx = -1
        right_side_sign = 1.0_8
    else
        call parcomm_global%abort("sbp_factory_mod, unknown sbp operator name: "// sbp_operator_name)
    end if

    nw_edge = size(sbp_op%W_edge_l,1)
    n_edge  = size(sbp_op%W_edge_l,2)
    sbp_op%W_edge_r = right_side_sign * sbp_op%W_edge_l(nw_edge:1:-1,n_edge:1:-1)
    sbp_op%edge_first_shift_r = -sbp_op%edge_last_l(n_edge:1:-1)+1

end function create_sbp_operator

end module sbp_factory_mod
