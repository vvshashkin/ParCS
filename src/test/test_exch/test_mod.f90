module test_mod

    use grid_field_mod,         only : grid_field_t
    use grid_field_factory_mod, only : create_grid_field

implicit none

contains

! subroutine test_cross_halo_vec_exchange()
!
!     use mpi
!
!     use exchange_vec_abstract_mod, only : exchange_vec_t
!     use exchange_vec_halo_mod,     only : exchange_vec_2D_halo_t
!     use partition_mod,             only : partition_t
!     use exchange_factory_mod,      only : create_Agrid_vec_halo_exchange
!     use topology_mod,              only : transform_index
!
!     class(exchange_vec_t), allocatable :: exch_vec_halo
!     type(partition_t)                  :: partition
!     type(grid_field_t) :: u, v
!
!     integer(kind=4)                    :: nh=100, nz=10, halo_width=50
!     integer(kind=4)                    :: myid, np, ierr, code
!
!     integer(kind=4) :: ts, te
!     integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum
!
!     integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
!
!     integer(kind=4)  :: i_out, j_out, i_step, j_step
!     character(len=1) :: first_dim_index
!
!     real(kind=8) :: u_remote, v_remote
!
!     logical lpass, glpass
!
!     call MPI_comm_rank(mpi_comm_world , myid, ierr)
!     call MPI_comm_size(mpi_comm_world , Np  , ierr)
!
!     if (myid==0) print*, 'Running vector full_halo_exchange test!'
!
!     call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')
!
!     !find start and end index of tiles belonging to the current proccesor
!     ts = partition%ts
!     te = partition%te
!
!     call create_grid_field(u, halo_width, 0, partition)
!     call create_grid_field(v, halo_width, 0, partition)
!
!     do t = ts, te
!         do k = partition%tile(t)%ks, partition%tile(t)%ke
!             do j = partition%tile(t)%js, partition%tile(t)%je
!                 do i = partition%tile(t)%is, partition%tile(t)%ie
!                     u%tile(t)%p(i,j,k) =  (partition%tile(t)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
!                     v%tile(t)%p(i,j,k) =  (partition%tile(t)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k + partition%num_panels*nh*nh*nz
!                 end do
!             end do
!         end do
!     end do
!
!     !Init exchange
!
!     exch_vec_halo = create_Agrid_vec_halo_exchange(partition, halo_width, 'full', myid, np)
!
!     !Perform exchange
!     call exch_vec_halo%do(u, v)
!
!     err_sum = 0
!     select type (exch_vec_halo)
!     class is (exchange_vec_2D_halo_t)
!     do ind = 1, exch_vec_halo%recv_number
!
!         local_tile_ind           = exch_vec_halo%recv_to_tile_ind(ind)
!         remote_tile_ind          = 1 + (exch_vec_halo%recv_tag(ind)-local_tile_ind)/(partition%num_panels*partition%num_tiles)
!
!         remote_tile_panel_number = partition%tile(remote_tile_ind)%panel_number
!         local_tile_panel_number  = partition%tile(local_tile_ind )%panel_number
!
!         do k = exch_vec_halo%recv_ks(ind), exch_vec_halo%recv_ke(ind)
!             do j = exch_vec_halo%recv_js(ind), exch_vec_halo%recv_je(ind)
!                 do i = exch_vec_halo%recv_is(ind), exch_vec_halo%recv_ie(ind)
!                     if (local_tile_panel_number == remote_tile_panel_number) then
!                         err_sum = err_sum + abs(int(u%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)))
!                     else
!                         call transform_index(remote_tile_panel_number, local_tile_panel_number, nh, i, j, i_out, j_out, i_step, j_step, first_dim_index)
!                         if (first_dim_index == "i") then
!                             u_remote = i_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
!                             v_remote = j_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k + partition%num_panels*nh*nh*nz)
!                             err_sum = err_sum + abs(int(u%tile(local_tile_ind)%p(i,j,k) - u_remote)) + abs(int(v%tile(local_tile_ind)%p(i,j,k) - v_remote))
!                         else
!                             v_remote = i_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
!                             u_remote = j_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k + partition%num_panels*nh*nh*nz)
!                             err_sum = err_sum + abs(int(u%tile(local_tile_ind)%p(i,j,k) - u_remote)) + abs(int(v%tile(local_tile_ind)%p(i,j,k) - v_remote))
!                         end if
!                     end if
!                 end do
!             end do
!         end do
!     end do
!     class default
!         call avost('Wrong type of vec halo exch object!')
!     end select
!     call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
!
!     if (gl_err_sum==0) then
!         if (myid==0) print*, 'Test passed!'
!     else
!         if (myid==0) print*, 'Test not passed! Error! Abort!'
!         call mpi_abort(mpi_comm_world, code, ierr)
!         stop
!     end if
!
!     contains
!
! end subroutine test_cross_halo_vec_exchange

subroutine test_A_halo_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_mod,      only : exchange_2D_halo_t
    use partition_mod,          only : partition_t
    use exchange_factory_mod,   only : create_symm_halo_exchange_A
    use topology_mod,           only : transform_index, find_basis_orientation

    type(domain_t) :: domain

    class(exchange_t), allocatable     :: exch_halo
    type(grid_field_t)                 :: f, u, v
    integer(kind=4)                    :: nh=100, nz=10, halo_width=50
    integer(kind=4)                    :: ierr, code

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote, f_remote

    call create_domain(domain, "cube", 'A', nh, nz)

    call domain%parcomm%print('Running cross_halo_exchange test!')

    call create_grid_field(f, halo_width, 0, domain%mesh_p)

    do t = domain%mesh_p%ts, domain%mesh_p%te
        do k = domain%mesh_p%tile(t)%ks, domain%mesh_p%tile(t)%ke
            do j = domain%mesh_p%tile(t)%js, domain%mesh_p%tile(t)%je
                do i = domain%mesh_p%tile(t)%is, domain%mesh_p%tile(t)%ie
                    f%tile(t)%p(i,j,k) = (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
                end do
            end do
        end do
    end do

    !Init exchange
    exch_halo = create_symm_halo_exchange_A(domain%partition, domain%parcomm, halo_width, 'full')

    !Perform exchange
    call exch_halo%do(f, domain%parcomm)

    err_sum = 0
    select type (exch_halo)
    class is (exchange_2D_halo_t)
    do ind = 1, exch_halo%recv_number

        local_tile_ind           = exch_halo%recv_to_tile_ind(ind)
        remote_tile_ind          = 1 + (exch_halo%recv_tag(ind)-local_tile_ind)/(domain%partition%num_panels*domain%partition%num_tiles)

        remote_tile_panel_number = domain%partition%panel_map(remote_tile_ind)
        local_tile_panel_number  = domain%partition%panel_map(local_tile_ind)

        do k = exch_halo%recv_tile(ind)%ks, exch_halo%recv_tile(ind)%ke
            do j = exch_halo%recv_tile(ind)%js, exch_halo%recv_tile(ind)%je
                do i = exch_halo%recv_tile(ind)%is, exch_halo%recv_tile(ind)%ie
                    call transform_index(local_tile_panel_number , i    , j    ,       &
                                         remote_tile_panel_number, i_out, j_out, nh, nh)
                    f_remote = ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
                    err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - f_remote))
                end do
            end do
        end do
    end do
    class default
        call avost('Wrong type of halo exch object!')
    end select
    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%print('Test failed! Error! Abort!')
        call mpi_abort(mpi_comm_world, code, ierr)
        stop
    end if

    call create_grid_field(u, halo_width, 0, domain%mesh_u)
    call create_grid_field(v, halo_width, 0, domain%mesh_v)


    do t = domain%mesh_u%ts, domain%mesh_u%te
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
                end do
            end do
        end do
    end do
    do t = domain%mesh_v%ts, domain%mesh_v%te
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k + domain%partition%num_panels*nh*nh*nz
                end do
            end do
        end do
    end do

    call exch_halo%do_vec(u, v, domain%parcomm)

        err_sum = 0
        select type (exch_halo)
        class is (exchange_2D_halo_t)
        do ind = 1, exch_halo%recv_number

            local_tile_ind  = exch_halo%recv_to_tile_ind(ind)
            remote_tile_ind = 1 + (exch_halo%recv_tag(ind)-local_tile_ind)/(domain%partition%num_panels*domain%partition%num_tiles)

            remote_tile_panel_number = domain%partition%panel_map(remote_tile_ind)
            local_tile_panel_number  = domain%partition%panel_map(local_tile_ind )

            do k = exch_halo%recv_tile(ind)%ks, exch_halo%recv_tile(ind)%ke
                do j = exch_halo%recv_tile(ind)%js, exch_halo%recv_tile(ind)%je
                    do i = exch_halo%recv_tile(ind)%is, exch_halo%recv_tile(ind)%ie
                        call transform_index(local_tile_panel_number , i    , j    ,       &
                                             remote_tile_panel_number, i_out, j_out, nh, nh)
                        call find_basis_orientation(local_tile_panel_number, remote_tile_panel_number, i_step, j_step, first_dim_index)
                        if (first_dim_index == "i") then
                            u_remote = i_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
                            v_remote = j_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k + domain%partition%num_panels*nh*nh*nz)
                        else
                            v_remote = i_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
                            u_remote = j_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k + domain%partition%num_panels*nh*nh*nz)
                        end if
                        err_sum = err_sum + abs(int(u%tile(local_tile_ind)%p(i,j,k) - u_remote)) + abs(int(v%tile(local_tile_ind)%p(i,j,k) - v_remote))
                    end do
                end do
            end do
        end do
        class default
            call avost('Wrong type of vec halo exch object!')
        end select
        call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

        if (gl_err_sum==0) then
            call domain%parcomm%print('Test passed!')
        else
            call domain%parcomm%abort('Test not passed! Error! Abort!')
        end if

end subroutine test_A_halo_exchange
subroutine test_halo_u_exchange()

    use mpi

    use domain_mod,             only : domain_t
    use domain_factory_mod,     only : create_domain
    use exchange_abstract_mod,  only : exchange_t
    use exchange_halo_C_mod,    only : exchange_2D_halo_C_t
    use partition_mod,          only : partition_t
    use exchange_factory_mod,   only : create_symm_halo_vec_exchange_U_points
    use topology_mod,           only : transform_index, find_basis_orientation

    type(domain_t) :: domain

    class(exchange_t), allocatable :: exch_halo
    type(grid_field_t)             :: f, u, v
    integer(kind=4)                :: nh=100, nz=10, halo_width=10
    integer(kind=4)                :: ierr

    integer(kind=4) :: ts, te
    integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum

    integer(kind=4)  :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
    integer(kind=4)  :: i_out, j_out, i_step, j_step
    character(len=1) :: first_dim_index

    real(kind=8) :: u_remote, v_remote

    call create_domain(domain, "cube", 'C', nh, nz)

    call domain%parcomm%print('Running u_halo_C_exchange test!')

    call create_grid_field(u, halo_width+1, 0, domain%mesh_u)
    call create_grid_field(v, halo_width+1, 0, domain%mesh_v)

    do t = domain%mesh_u%ts, domain%mesh_u%te
        do k = domain%mesh_u%tile(t)%ks, domain%mesh_u%tile(t)%ke
            do j = domain%mesh_u%tile(t)%js, domain%mesh_u%tile(t)%je
                do i = domain%mesh_u%tile(t)%is, domain%mesh_u%tile(t)%ie
                    u%tile(t)%p(i,j,k) =  (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
                end do
            end do
        end do
    end do
    do t = domain%mesh_v%ts, domain%mesh_v%te
        do k = domain%mesh_v%tile(t)%ks, domain%mesh_v%tile(t)%ke
            do j = domain%mesh_v%tile(t)%js, domain%mesh_v%tile(t)%je
                do i = domain%mesh_v%tile(t)%is, domain%mesh_v%tile(t)%ie
                    v%tile(t)%p(i,j,k) =  (domain%partition%panel_map(t)-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k + domain%partition%num_panels*nh*nh*nz
                end do
            end do
        end do
    end do

    !Init exchange

    exch_halo = create_symm_halo_vec_exchange_U_points(domain%partition, domain%parcomm, halo_width, 'full')

    !Perform exchange
    call exch_halo%do_vec(u, v, domain%parcomm)

    err_sum = 0
    select type (exch_halo)
    class is (exchange_2D_halo_C_t)
    do ind = 1, exch_halo%recv_number

        local_tile_ind           = exch_halo%recv_to_tile_ind(ind)
        remote_tile_ind          = 1 + (exch_halo%recv_tag(ind)-local_tile_ind)/(domain%partition%num_panels*domain%partition%num_tiles)

        remote_tile_panel_number = domain%partition%panel_map(remote_tile_ind)
        local_tile_panel_number  = domain%partition%panel_map(local_tile_ind )

        do k = exch_halo%recv_tile(ind)%ks, exch_halo%recv_tile(ind)%ke
            do j = exch_halo%recv_tile(ind)%js, exch_halo%recv_tile(ind)%je
                do i = exch_halo%recv_tile(ind)%is, exch_halo%recv_tile(ind)%ie
                    call find_basis_orientation(local_tile_panel_number, remote_tile_panel_number, i_step, j_step, first_dim_index)
                        if (first_dim_index == "i") then
                            call transform_index(local_tile_panel_number, i, j, remote_tile_panel_number, i_out, j_out, nh+1, nh)
                            u_remote = i_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k)
                        else
                            call transform_index(local_tile_panel_number, i, j, remote_tile_panel_number, i_out, j_out, nh, nh+1)
                            u_remote = j_step*((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j_out-1) + nz*(i_out-1) + k + domain%partition%num_panels*nh*nh*nz)
                        end if
                        err_sum = err_sum + abs(int(u%tile(local_tile_ind)%p(i,j,k) - u_remote))
                end do
            end do
        end do
    end do
    class default
        call domain%parcomm%abort('Wrong type of vec halo exch object!')
    end select

    call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, domain%parcomm%comm_w, ierr)

    if (gl_err_sum==0) then
        call domain%parcomm%print('Test passed!')
    else
        call domain%parcomm%abort('Test not passed! Error! Abort!')
    end if
!
end subroutine test_halo_u_exchange

! subroutine test_full_halo_exchange()
!
!     use mpi
!
!     use exchange_abstract_mod, only : exchange_t
!     use exchange_halo_mod,     only : exchange_2D_halo_t
!     use partition_mod,         only : partition_t
!     use exchange_factory_mod,  only : create_Agrid_halo_exchange
!
!     class(exchange_t), allocatable     :: exch_halo
!     type(partition_t)                  :: partition
!     type(grid_field_t) :: f
!
!     integer(kind=4)                    :: nh=100, nz=10, halo_width=50
!     integer(kind=4)                    :: myid, np, ierr, code
!
!     integer(kind=4) :: ts, te
!     integer(kind=4) :: ind, t, i, j, k, err_sum, gl_err_sum
!
!     integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
!
!     call MPI_comm_rank(mpi_comm_world , myid, ierr)
!     call MPI_comm_size(mpi_comm_world , Np  , ierr)
!
!     if (myid==0) print*, 'Running cross_halo_exchange test!'
!
!     call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')
!
!     !find start and end index of tiles belonging to the current proccesor
!     ts = partition%ts
!     te = partition%te
!
!     !Init arrays
!     call create_grid_field(f, halo_width, 0, partition)
!
!     do t = ts, te
!         do k = partition%tile(t)%ks, partition%tile(t)%ke
!             do j = partition%tile(t)%js, partition%tile(t)%je
!                 do i = partition%tile(t)%is, partition%tile(t)%ie
!                     f%tile(t)%p(i,j,k) =  (partition%tile(t)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
!                 end do
!             end do
!         end do
!     end do
!
!     !Init exchange
!
!     exch_halo = create_Agrid_halo_exchange(partition, halo_width, 'full', myid, np)
!
!     !Perform exchange
!     call exch_halo%do(f)
!
!     err_sum = 0
!     select type (exch_halo)
!     class is (exchange_2D_halo_t)
!     do ind = 1, exch_halo%recv_number
!
!         local_tile_ind           = exch_halo%recv_to_tile_ind(ind)
!         remote_tile_ind          = 1 + (exch_halo%recv_tag(ind)-local_tile_ind)/(6*partition%num_tiles)
!
!         remote_tile_panel_number = partition%tile(remote_tile_ind)%panel_number
!         local_tile_panel_number  = partition%tile(local_tile_ind )%panel_number
!
!         do k = exch_halo%recv_ks(ind), exch_halo%recv_ke(ind)
!             do j = exch_halo%recv_js(ind), exch_halo%recv_je(ind)
!                 do i = exch_halo%recv_is(ind), exch_halo%recv_ie(ind)
!
!                     if (local_tile_panel_number == remote_tile_panel_number) then
!                         err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k)))
!                     else
!                         if (exch_halo%send_j_step(ind)==1 .and. exch_halo%send_i_step(ind)==1 .and. exch_halo%first_dim_index(ind)=='i') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(modulo(i-1,nh)+1-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==1 .and. exch_halo%send_i_step(ind)==1 .and. exch_halo%first_dim_index(ind)=='j') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(i-1,nh)+1-1) + nz*(modulo(j-1,nh)+1-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==1 .and. exch_halo%send_i_step(ind)==-1 .and. exch_halo%first_dim_index(ind)=='i') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(modulo(j-1,nh)+1-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==1 .and. exch_halo%send_i_step(ind)==-1 .and. exch_halo%first_dim_index(ind)=='j') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(modulo(j-1,nh)+1-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==-1 .and. exch_halo%send_i_step(ind)==1 .and. exch_halo%first_dim_index(ind)=='i') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(modulo(i-1,nh)+1-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==-1 .and. exch_halo%send_i_step(ind)==1 .and. exch_halo%first_dim_index(ind)=='j') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(modulo(i-1,nh)+1-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==-1 .and. exch_halo%send_i_step(ind)==-1 .and. exch_halo%first_dim_index(ind)=='i') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*nh*(nh-modulo(j-1,nh)-1) + nz*(nh-modulo(i-1,nh)-1) + k)))
!                         else if (exch_halo%send_j_step(ind)==-1 .and. exch_halo%send_i_step(ind)==-1 .and. exch_halo%first_dim_index(ind)=='j') then
!                             err_sum = err_sum + abs(int(f%tile(local_tile_ind)%p(i,j,k) - ((remote_tile_panel_number-1)*nh*nh*nz + nz*(nh-modulo(j-1,nh)-1) + nz*nh*(nh-modulo(i-1,nh)-1) + k)))
!                         else
!                             print*, 'Error!!!', 'myid=', myid
!                             call mpi_abort(mpi_comm_world, code, ierr)
!                             stop
!                         end if
!                     end if
!                 end do
!             end do
!         end do
!     end do
!     class default
!         call avost('Wrong type of halo exch object!')
!     end select
!     call mpi_allreduce(err_sum, gl_err_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
!
!     if (gl_err_sum==0) then
!         if (myid==0) print*, 'Test passed!'
!     else
!         if (myid==0) print*, 'Test not passed! Error! Abort!'
!         call mpi_abort(mpi_comm_world, code, ierr)
!         stop
!     end if
!
! end subroutine test_full_halo_exchange

! subroutine test_gather_exchange()
!
! use mpi
!
! use exchange_abstract_mod, only : exchange_t
! use partition_mod,         only : partition_t
! use exchange_factory_mod,  only : create_gather_exchange
!
! class(exchange_t),allocatable      :: exch_gather
! type(partition_t)                  :: partition
! type(grid_field_t)                 :: f
!
! integer(kind=4)                    :: nh=100, nz=10, halo_width=10
! integer(kind=4)                    :: myid, np, ierr, code
! integer(kind=4)                    :: master_id = 0
!
! integer(kind=4) :: ts, te
! integer(kind=4) :: ind, i, j, k, t
!
! real(kind=8) :: err_sum
!
! integer(kind=4) :: local_tile_ind, remote_tile_ind, local_tile_panel_number, remote_tile_panel_number
!
! call MPI_comm_rank(mpi_comm_world , myid, ierr)
! call MPI_comm_size(mpi_comm_world , Np  , ierr)
!
!
! if (myid==0) print*, 'Running gather_exchange test!'
!
! call partition%init(nh, nz, max(1,Np/6), myid, Np, strategy = 'default')
!
! !find start and end index of tiles belonging to the current proccesor
! ts = partition%ts
! te = partition%te
!
! !Init arrays
!
! if (myid == master_id) then
!     allocate(f%tile(1:6*partition%num_tiles))
!     do i = 1, 6*partition%num_tiles
!         call f%tile(i)%init(partition%tile(i)%panel_number,             &
!                         partition%tile(i)%is, partition%tile(i)%ie, &
!                         partition%tile(i)%js, partition%tile(i)%je, &
!                         partition%tile(i)%ks, partition%tile(i)%ke, &
!                         halo_width, halo_width, 0)
!         f%tile(i)%p = huge(1.0_8)
!     end do
! else
!     call create_grid_field(f, halo_width, 0, partition)
! end if
!
! do t = ts, te
!     do k = partition%tile(t)%ks, partition%tile(t)%ke
!         do j = partition%tile(t)%js, partition%tile(t)%je
!             do i = partition%tile(t)%is, partition%tile(t)%ie
!                 f%tile(t)%p(i,j,k) =  (partition%tile(t)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k
!             end do
!         end do
!     end do
! end do
!
! !Init exchange
! exch_gather = create_gather_exchange(partition, master_id, myid, np)
!
! !Perform exchange
! call exch_gather%do(f)
!
! if (myid == master_id) then
!
!     err_sum = 0
!
!     do ind = 1, 6*partition%num_tiles
!         do k = partition%tile(ind)%ks, partition%tile(ind)%ke
!             do j = partition%tile(ind)%js, partition%tile(ind)%je
!                 do i = partition%tile(ind)%is, partition%tile(ind)%ie
!                     err_sum = err_sum + abs(f%tile(ind)%p(i,j,k) - ((partition%tile(ind)%panel_number-1)*nh*nh*nz + nz*nh*(j-1) + nz*(i-1) + k))
!                 end do
!             end do
!         end do
!     end do
!
!     if (int(err_sum)==0) then
!         print*, 'Test passed!'
!     else
!         print*, 'Test not passed! Error! Abort!'
!         call mpi_abort(mpi_comm_world, code, ierr)
!         stop
!     end if
!
! end if
!
! end subroutine test_gather_exchange


end module test_mod
