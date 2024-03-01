program mpi_potts_simulation
  use, intrinsic :: iso_fortran_env
  use mpi
  use potts_m
  use kahan_summation_m
  use mpi_variance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000
  integer(int32), parameter :: nsample = 25
  integer(int64), parameter :: nx = 1001
  integer(int64), parameter :: ny = 1000
  integer(int32), parameter :: state = 5
  real(real64), parameter :: kbt = 2.26918531421302_real64
  real(real64), parameter :: n_inv_r64 = 1 / real(nx * ny, real64)
  type(potts) :: system
  type(variance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, j
  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  call system%init(nx, ny, kbt, state)
  if (myrank == 0) then
     write(output_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
     write(output_unit, '(3(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny(), " state: ", system%max_state()
     write(output_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(output_unit, '(a, g0)' ) "# 温度: ", system%kbt()
     write(output_unit, '(a)' ) "# method: Metropolis"
     write(output_unit, '(a, i0)' ) "# the number of processors: ", num_proc

     write(error_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
     write(error_unit, '(2(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny(), " state: ", system%max_state()
     write(error_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(error_unit, '(a, g0)' ) "# 温度: ", system%kbt()
     write(error_unit, '(a)' ) "# method: Metropolis"
     write(error_unit, '(a, i0)' ) "# the number of processors: ", num_proc
  end if
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(a, i0)') "sample: ", j
     call system%set_ising_allup()
     do i = 1, mcs
        call system%update()
        associate(m => system%calc_magne_summ())
          call order_parameter(i)%add_data(m * n_inv_r64)
        end associate
     end do
  end do
  block
    type(variance_kahan) :: all_order_params(mcs)
    call vk_mpi_multi_gather(mcs, order_parameter, all_order_params, 0, myrank, num_proc, ierr)
    if (myrank == 0) then
       write(output_unit, '(a)') "# Nsize, Nsample, mcs, <m>, <m^2>, χ'"
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))') system%nall(), all_order_params(i)%num_sample(), i, &
               & all_order_params(i)%mean(), all_order_params(i)%square_mean(), system%nall() * all_order_params(i)%var()
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_potts_simulation
