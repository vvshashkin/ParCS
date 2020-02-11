module exchange_factory_mod

use partition_mod,        only : partition_t
use tile_mod,             only : tile_t

implicit none

contains

function create_2d_halo_exchange(partition, halo_width, halo_type, myid, np ) result(exchange)

    use exchange_halo_mod,    only : exchange_2D_halo_t

    type(partition_t), target, intent(in)    :: partition

    integer(kind=4), intent(in) :: halo_width
    character(*),    intent(in) :: halo_type
    integer(kind=4), intent(in) :: myid, np

    type(exchange_2D_halo_t), allocatable :: exchange

    type(tile_t), pointer :: local_tile, remote_tile

    integer(kind=4),  dimension((6*partition%num_tiles)**2) :: recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke, &
                                                               send_is, send_ie, send_js, send_je, send_ks, send_ke
    integer(kind=4),  dimension((6*partition%num_tiles)**2) :: send_to_proc_id, recv_from_proc_id   , &
                                                               recv_to_tile_ind , &
                                                               send_from_tile_ind, &
                                                               send_tag, recv_tag
    integer(kind=4),  dimension((6*partition%num_tiles)**2) :: send_pts_num, recv_pts_num

    integer(kind=4),  dimension((6*partition%num_tiles)**2) :: send_i_step, send_j_step
    character(len=1), dimension((6*partition%num_tiles)**2) :: first_dim_index


    integer(kind=4) :: local_ind, remote_ind, exch_num, ind,ts ,te
    logical :: is_intersection

    if (halo_width > partition%Nh) then
        write(*,*) 'Error! Halo is too wide! Abort!'
        stop
    end if

    ts = findloc(partition%proc_map, myid, dim=1)
    te = findloc(partition%proc_map, myid, back = .true., dim=1)

    exch_num = 0

    do local_ind = 1, 6*partition%num_tiles

        if (partition%proc_map(local_ind) /= myid) cycle

        local_tile => partition%tile(local_ind)

        do remote_ind = 1, 6*partition%num_tiles

            remote_tile => partition%tile(remote_ind)

            select case (halo_type)

            case('cross')

            call find_cross_halo_intersection(local_tile, remote_tile, halo_width, partition%Nh, is_intersection , &
                                              recv_is(exch_num + 1), recv_ie(exch_num + 1)                            , &
                                              recv_js(exch_num + 1), recv_je(exch_num + 1)                            , &
                                              send_is(exch_num + 1), send_ie(exch_num + 1)                            , &
                                              send_js(exch_num + 1), send_je(exch_num + 1)                            , &
                                              send_i_step(exch_num+1), send_j_step(exch_num+1), first_dim_index(exch_num+1) )

            case('full')
                call find_full_halo_intersection(local_tile, remote_tile, halo_width, partition%Nh, is_intersection , &
                                                  recv_is(exch_num + 1), recv_ie(exch_num + 1)                            , &
                                                  recv_js(exch_num + 1), recv_je(exch_num + 1)                            , &
                                                  send_is(exch_num + 1), send_ie(exch_num + 1)                            , &
                                                  send_js(exch_num + 1), send_je(exch_num + 1)                            , &
                                                  send_i_step(exch_num+1), send_j_step(exch_num+1), first_dim_index(exch_num+1) )
            case default
                call avost("Wrong halo_type in create_2d_halo_exchange function! Stop!")
            end select

            if (is_intersection) then
                exch_num = exch_num + 1

                recv_to_tile_ind  (exch_num) = local_ind

                send_from_tile_ind(exch_num) = local_ind

                send_ks(exch_num) = local_tile%ks
                send_ke(exch_num) = local_tile%ke

                recv_ks(exch_num) = remote_tile%ks
                recv_ke(exch_num) = remote_tile%ke

                send_pts_num(exch_num) = (send_ie(exch_num) - send_is(exch_num) + 1)* &
                                         (send_je(exch_num) - send_js(exch_num) + 1)* &
                                         (send_ke(exch_num) - send_ks(exch_num) + 1)

                recv_pts_num(exch_num) = (recv_ie(exch_num) - recv_is(exch_num) + 1)* &
                                         (recv_je(exch_num) - recv_js(exch_num) + 1)* &
                                         (recv_ke(exch_num) - recv_ks(exch_num) + 1)

                send_to_proc_id  (exch_num) = partition%proc_map(remote_ind)
                recv_from_proc_id(exch_num) = partition%proc_map(remote_ind)

                send_tag(exch_num) = 6*partition%num_tiles*(local_ind-1) + remote_ind
                recv_tag(exch_num) = 6*partition%num_tiles*(remote_ind-1) + local_ind

            end if

        end do

    end do

    allocate(exchange, source = exchange_2D_halo_t(              &
            send_is            = send_is(1:exch_num),          &
            send_ie            = send_ie(1:exch_num),          &
            send_js            = send_js(1:exch_num),          &
            send_je            = send_je(1:exch_num),          &
            send_ks            = send_ks(1:exch_num),          &
            send_ke            = send_ke(1:exch_num),          &
            recv_is            = recv_is(1:exch_num),          &
            recv_ie            = recv_ie(1:exch_num),          &
            recv_js            = recv_js(1:exch_num),          &
            recv_je            = recv_je(1:exch_num),          &
            recv_ks            = recv_ks(1:exch_num),          &
            recv_ke            = recv_ke(1:exch_num),          &
            send_i_step        = send_i_step(1:exch_num),      &
            send_j_step        = send_j_step(1:exch_num),      &
            first_dim_index    = first_dim_index(1:exch_num),  &
            recv_to_tile_ind   = recv_to_tile_ind(1:exch_num),    &
            send_from_tile_ind = send_from_tile_ind(1:exch_num),    &
            send_to_proc_id    = send_to_proc_id(1:exch_num),     &
            recv_from_proc_id  = recv_from_proc_id(1:exch_num),     &
            send_points_num    = send_pts_num(1:exch_num),     &
            recv_points_num    = recv_pts_num(1:exch_num),     &
            send_tag           = send_tag(1:exch_num),         &
            recv_tag           = recv_tag(1:exch_num),         &
            send_number        = exch_num                ,     &
            recv_number        = exch_num, &
            ts                 = ts, &
            te                 = te) )

    ! call exchange%profile%check()
    !
     allocate(exchange%send_buff(exch_num)   , exchange%recv_buff(exch_num)   )
     allocate(exchange%mpi_send_req(exch_num), exchange%mpi_recv_req(exch_num))
    do ind = 1, exch_num
        call exchange%send_buff(ind)%init(send_pts_num(ind))
        call exchange%recv_buff(ind)%init(recv_pts_num(ind))
    end do


end function create_2d_halo_exchange

subroutine find_cross_halo_intersection(l_tile, r_tile, halo_width, npoints, is_intersection, &
                                        recv_is, recv_ie, recv_js, recv_je,                   &
                                        send_is, send_ie, send_js, send_je,                   &
                                        i_step, j_step, fisrt_dim)

    use topology_mod, only : transform_index

    type(tile_t),     intent(in)  :: l_tile, r_tile
    integer(kind=4),  intent(in)  :: halo_width, npoints
    logical,          intent(out) :: is_intersection
    integer(kind=4),  intent(out) :: recv_is, recv_ie, recv_js, recv_je, &
                                    send_is, send_ie, send_js, send_je
    integer(kind=4),  intent(out) :: i_step, j_step
    character(len=1), intent(out) :: fisrt_dim

    integer(kind=4)  :: r_tile_is, r_tile_ie, r_tile_js, r_tile_je
    integer(kind=4)  :: is_r, ie_r, js_r, je_r
    integer(kind=4)  :: is_l, ie_l, js_l, je_l
    integer(kind=4)  :: vertex_i(4), vertex_j(4)
    is_intersection = .true.

    if (l_tile%panel_number == r_tile%panel_number) then
        r_tile_is = r_tile%is; r_tile_ie = r_tile%ie;
        r_tile_js = r_tile%js; r_tile_je = r_tile%je;
        i_step = 1
        j_step = 1
        fisrt_dim = 'i'
    else
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%js, vertex_i(1), vertex_j(1))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%je, vertex_i(2), vertex_j(2))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%js, vertex_i(3), vertex_j(3))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%je, vertex_i(4), vertex_j(4), i_step, j_step, fisrt_dim)

        r_tile_is = minval(vertex_i); r_tile_ie = maxval(vertex_i)
        r_tile_js = minval(vertex_j); r_tile_je = maxval(vertex_j)

    end if

    !left recv halo
    is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
    js_l = l_tile%js           ; je_l = l_tile%je

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !left send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_ie+1; ie_r = r_tile_ie+halo_width
        js_r = r_tile_js  ; je_r = r_tile_je

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !right recv halo
    is_l = l_tile%ie+1; ie_l = l_tile%ie+halo_width
    js_l = l_tile%js  ; je_l = l_tile%je

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !right send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is - halo_width; ie_r = r_tile_is - 1
        js_r = r_tile_js             ; je_r = r_tile_je

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !bottom recv halo
    is_l = l_tile%is           ; ie_l = l_tile%ie
    js_l = l_tile%js-halo_width; je_l = l_tile%js-1

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !bottom send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is    ; ie_r = r_tile_ie
        js_r = r_tile_je + 1; je_r = r_tile_je + halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !top recv halo
    is_l = l_tile%is  ; ie_l = l_tile%ie
    js_l = l_tile%je+1; je_l = l_tile%je+halo_width

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !top send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is             ; ie_r = r_tile_ie
        js_r = r_tile_js - halo_width; je_r = r_tile_js - 1

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    is_intersection = .false.

end subroutine find_cross_halo_intersection

subroutine find_full_halo_intersection(l_tile, r_tile, halo_width, npoints, is_intersection, &
                                        recv_is, recv_ie, recv_js, recv_je,                   &
                                        send_is, send_ie, send_js, send_je,                   &
                                        i_step, j_step, fisrt_dim)

    use topology_mod, only : transform_index

    type(tile_t),     intent(in)  :: l_tile, r_tile
    integer(kind=4),  intent(in)  :: halo_width, npoints
    logical,          intent(out) :: is_intersection
    integer(kind=4),  intent(out) :: recv_is, recv_ie, recv_js, recv_je, &
                                    send_is, send_ie, send_js, send_je
    integer(kind=4),  intent(out) :: i_step, j_step
    character(len=1), intent(out) :: fisrt_dim

    integer(kind=4)  :: r_tile_is, r_tile_ie, r_tile_js, r_tile_je
    integer(kind=4)  :: is_r, ie_r, js_r, je_r
    integer(kind=4)  :: is_l, ie_l, js_l, je_l
    integer(kind=4)  :: vertex_i(4), vertex_j(4)
    is_intersection = .true.

    if (l_tile%panel_number == r_tile%panel_number) then
        r_tile_is = r_tile%is; r_tile_ie = r_tile%ie;
        r_tile_js = r_tile%js; r_tile_je = r_tile%je;
        i_step = 1
        j_step = 1
        fisrt_dim = 'i'
    else
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%js, vertex_i(1), vertex_j(1))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%is, r_tile%je, vertex_i(2), vertex_j(2))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%js, vertex_i(3), vertex_j(3))
        call transform_index(l_tile%panel_number, r_tile%panel_number, npoints, r_tile%ie, r_tile%je, vertex_i(4), vertex_j(4), i_step, j_step, fisrt_dim)

        r_tile_is = minval(vertex_i); r_tile_ie = maxval(vertex_i)
        r_tile_js = minval(vertex_j); r_tile_je = maxval(vertex_j)

    end if

    !left recv halo
    is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
    js_l = l_tile%js           ; je_l = l_tile%je

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !left send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
        js_l = l_tile%js-halo_width; je_l = l_tile%je+halo_width

        recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
        recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_ie+1           ; ie_r = r_tile_ie+halo_width
        js_r = r_tile_js-halo_width  ; je_r = r_tile_je+halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !right recv halo
    is_l = l_tile%ie+1; ie_l = l_tile%ie+halo_width
    js_l = l_tile%js  ; je_l = l_tile%je

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !right send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%ie+1          ; ie_l = l_tile%ie+halo_width
        js_l = l_tile%js-halo_width ; je_l = l_tile%je+halo_width

        recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
        recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is - halo_width; ie_r = r_tile_is - 1
        js_r = r_tile_js - halo_width; je_r = r_tile_je + halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !bottom recv halo
    is_l = l_tile%is           ; ie_l = l_tile%ie
    js_l = l_tile%js-halo_width; je_l = l_tile%js-1

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !bottom send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is-halo_width; ie_l = l_tile%ie+halo_width
        js_l = l_tile%js-halo_width; je_l = l_tile%js-1

        recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
        recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is-halo_width; ie_r = r_tile_ie + halo_width
        js_r = r_tile_je + 1       ; je_r = r_tile_je + halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !top recv halo
    is_l = l_tile%is  ; ie_l = l_tile%ie
    js_l = l_tile%je+1; je_l = l_tile%je+halo_width

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !top send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is-halo_width ; ie_l = l_tile%ie+halo_width
        js_l = l_tile%je+1          ; je_l = l_tile%je+halo_width

        recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
        recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is - halo_width; ie_r = r_tile_ie + halo_width
        js_r = r_tile_js - halo_width; je_r = r_tile_js - 1

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !bottom-left recv halo
    is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
    js_l = l_tile%js-halo_width; je_l = l_tile%js-1

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !bottom-left send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_ie+1  ; ie_r = r_tile_ie+halo_width
        js_r = r_tile_je+1  ; je_r = r_tile_je+halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !bottom-right recv halo
    is_l = l_tile%ie+1         ; ie_l = l_tile%ie+halo_width
    js_l = l_tile%js-halo_width; je_l = l_tile%js-1

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !bottom-right send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is - halo_width; ie_r = r_tile_is - 1
        js_r = r_tile_je+1           ; je_r = r_tile_je+halo_width

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !top-left recv halo
    is_l = l_tile%is-halo_width; ie_l = l_tile%is-1
    js_l = l_tile%je+1         ; je_l = l_tile%je+halo_width

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !top-left send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_ie+1           ; ie_r = r_tile_ie+halo_width
        js_r = r_tile_js-halo_width  ; je_r = r_tile_js-1

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    !top-right recv halo
    is_l = l_tile%ie+1         ; ie_l = l_tile%ie+halo_width
    js_l = l_tile%je+1         ; je_l = l_tile%je+halo_width

    recv_is = max(is_l, r_tile_is); recv_ie = min(ie_l, r_tile_ie)
    recv_js = max(js_l, r_tile_js); recv_je = min(je_l, r_tile_je)

    !top-right send halo
    if ( (recv_is<=recv_ie) .and. (recv_js<=recv_je) ) then

        is_l = l_tile%is; ie_l = l_tile%ie
        js_l = l_tile%js; je_l = l_tile%je

        is_r = r_tile_is - halo_width ; ie_r = r_tile_is - 1
        js_r = r_tile_js - halo_width ; je_r = r_tile_js-1

        send_is = max(is_l, is_r); send_ie = min(ie_l, ie_r)
        send_js = max(js_l, js_r); send_je = min(je_l, je_r)

        return
    end if

    is_intersection = .false.

end subroutine find_full_halo_intersection

function create_gather_exchange(partition, master_id, myid, np) result(exchange)

    use exchange_gather_mod,    only : exchange_gather_t

    type(partition_t), target, intent(in)    :: partition
    integer(kind=4),           intent(in)    :: master_id
    integer(kind=4),           intent(in)    :: myid, np

    type(exchange_gather_t), allocatable     :: exchange

    integer(kind=4),  dimension((6*partition%num_tiles)) :: recv_is, recv_ie, recv_js, recv_je, recv_ks, recv_ke, &
                                                                send_is, send_ie, send_js, send_je, send_ks, send_ke
    integer(kind=4),  dimension((6*partition%num_tiles)) :: recv_from_proc_id , recv_to_tile_ind , &
                                                                send_from_tile_ind,  send_tag, recv_tag
    integer(kind=4),  dimension((6*partition%num_tiles)) :: send_pts_num, recv_pts_num
    integer(kind=4) :: ind, send_exch_num, recv_exch_num, ts, te

    ts = findloc(partition%proc_map, myid, dim=1)
    te = findloc(partition%proc_map, myid, back = .true., dim=1)

    recv_exch_num = 0
    send_exch_num = 0

    if (myid == master_id) then


        do ind = 1, 6*partition%num_tiles

            if (partition%proc_map(ind) == myid) cycle

            recv_exch_num = recv_exch_num + 1

            recv_is(recv_exch_num) = partition%tile(ind)%is
            recv_ie(recv_exch_num) = partition%tile(ind)%ie

            recv_js(recv_exch_num) = partition%tile(ind)%js
            recv_je(recv_exch_num) = partition%tile(ind)%je

            recv_ks(recv_exch_num) = partition%tile(ind)%ks
            recv_ke(recv_exch_num) = partition%tile(ind)%ke

            recv_pts_num(recv_exch_num) = (recv_ie(recv_exch_num) - recv_is(recv_exch_num) + 1)* &
                                          (recv_je(recv_exch_num) - recv_js(recv_exch_num) + 1)* &
                                          (recv_ke(recv_exch_num) - recv_ks(recv_exch_num) + 1)

            recv_to_tile_ind  (recv_exch_num) = ind

            recv_from_proc_id(recv_exch_num)  = partition%proc_map(ind)

            recv_tag(recv_exch_num) = ind

        end do

    else

        do ind = 1, 6*partition%num_tiles

            if (partition%proc_map(ind) /= myid) cycle

            send_exch_num = send_exch_num + 1

            send_is(send_exch_num) = partition%tile(ind)%is
            send_ie(send_exch_num) = partition%tile(ind)%ie

            send_js(send_exch_num) = partition%tile(ind)%js
            send_je(send_exch_num) = partition%tile(ind)%je

            send_ks(send_exch_num) = partition%tile(ind)%ks
            send_ke(send_exch_num) = partition%tile(ind)%ke

            send_pts_num(send_exch_num) = (send_ie(send_exch_num) - send_is(send_exch_num) + 1)* &
                                          (send_je(send_exch_num) - send_js(send_exch_num) + 1)* &
                                          (send_ke(send_exch_num) - send_ks(send_exch_num) + 1)

            send_from_tile_ind(send_exch_num) = ind

            send_tag(send_exch_num) = ind

        end do
    end if

    allocate(exchange,         source = exchange_gather_t(            &
            send_is            = send_is(1:send_exch_num),            &
            send_ie            = send_ie(1:send_exch_num),            &
            send_js            = send_js(1:send_exch_num),            &
            send_je            = send_je(1:send_exch_num),            &
            send_ks            = send_ks(1:send_exch_num),            &
            send_ke            = send_ke(1:send_exch_num),            &
            recv_is            = recv_is(1:recv_exch_num),            &
            recv_ie            = recv_ie(1:recv_exch_num),            &
            recv_js            = recv_js(1:recv_exch_num),            &
            recv_je            = recv_je(1:recv_exch_num),            &
            recv_ks            = recv_ks(1:recv_exch_num),            &
            recv_ke            = recv_ke(1:recv_exch_num),            &
            recv_to_tile_ind   = recv_to_tile_ind(1:recv_exch_num),   &
            send_from_tile_ind = send_from_tile_ind(1:send_exch_num), &
            recv_from_proc_id  = recv_from_proc_id(1:recv_exch_num),  &
            send_points_num    = send_pts_num(1:send_exch_num),       &
            recv_points_num    = recv_pts_num(1:recv_exch_num),       &
            send_tag           = send_tag(1:send_exch_num),           &
            recv_tag           = recv_tag(1:recv_exch_num),           &
            send_number        = send_exch_num,                       &
            recv_number        = recv_exch_num,                       &
            master_id = master_id                                 ,   &
            ts = ts,                                                  &
            te = te ))

    allocate(exchange%send_buff(1:send_exch_num)   , exchange%recv_buff(1:recv_exch_num)   )
    allocate(exchange%mpi_send_req(1:send_exch_num), exchange%mpi_recv_req(1:recv_exch_num))
    do ind = 1, recv_exch_num
        call exchange%recv_buff(ind)%init(recv_pts_num(ind))
    end do
    do ind = 1, send_exch_num
        call exchange%send_buff(ind)%init(send_pts_num(ind))
    end do

end function create_gather_exchange

end module exchange_factory_mod
