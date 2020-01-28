module outputer_factory_mod

implicit none

contains

function create_master_process_outputer(master_id, gather_exch, write_type) result(master_process_outputer)

    use master_process_outputer_mod, only : master_process_outputer_t
    use exchange_mod,                only : exchange_t

    integer(kind=4),                 intent(in) :: master_id
    type(exchange_t),                intent(in) :: gather_exch
    character(*),                    intent(in) :: write_type ! 'bin' or 'txt' output

    type(master_process_outputer_t)             :: master_process_outputer

    master_process_outputer%master_id   = master_id
    master_process_outputer%gather_exch = gather_exch
    master_process_outputer%write_type  = write_type

end function create_master_process_outputer

end module outputer_factory_mod
