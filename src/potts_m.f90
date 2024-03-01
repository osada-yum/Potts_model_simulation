module potts_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: potts
  integer(int32), parameter :: potts_max_state_limit = 50 ! 実装的に 50^5 が限界.
  type :: potts
     private
     integer(int32) :: max_state_
     integer(int64) :: nx_, ny_, nall_
     real(real64) :: beta_
     integer(int32), allocatable :: spins_(:)
     real(real64), allocatable :: ws_(:, :, :, :, :, :)
   contains
     !> initializer.
     procedure, pass :: init => init_potts
     !> setter.
     procedure, pass :: set_ising_allup => set_ising_allup_potts
     procedure, pass :: set_ising_mixed_phase => set_ising_mixed_phase_potts
     procedure, pass :: set_ising_random => set_ising_random_potts
     procedure, pass :: set_kbt => set_kbt_potts
     procedure, pass :: set_beta => set_beta_potts
     !> updater.
     procedure, pass :: update => update_potts
     procedure, pass, private :: update_onesite => update_onesite_potts
     procedure, pass, private :: update_norishiro => update_norishiro_potts
     !> calculator.
     procedure, pass :: calc_energy_summ => calc_total_energy_potts
     procedure, pass :: calc_magne_summ => calc_total_magne_potts
     !> getter.
     procedure, pass :: nx => nx_potts
     procedure, pass :: ny => ny_potts
     procedure, pass :: nall => nall_potts
     procedure, pass :: max_state => max_state_potts
     procedure, pass :: kbt => kbt_potts
     procedure, pass :: beta => beta_potts
     procedure, pass, private :: norishiro_begin => norishiro_begin_potts
     procedure, pass, private :: norishiro_end => norishiro_end_potts
  end type potts
contains
  !> init_potts: Initialize potts once.
  pure subroutine init_potts(this, nx, ny, kbt, state)
    class(potts), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    integer(int32), intent(in) :: state
    if (allocated(this%spins_)) return
    if (is_even(nx) .or. (.not. is_even(ny))) then
       error stop "The parity of size must be (x, y) == (odd, even)."
    end if
    if (state > potts_max_state_limit) then
       error stop "The state of Potts must be smaller than 50 for this object."
    end if
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    this%max_state_ = state
    allocate(this%spins_(this%norishiro_begin() : this%norishiro_end()))
    call this%set_ising_allup()
    allocate(this%ws_(state, state, state, state, state, state))
    call this%set_kbt(kbt)
  contains
    pure logical function is_even(v) result(res)
      integer(int64), intent(in) :: v
      res = iand(v, b'1') == 0_int64
    end function is_even
  end subroutine init_potts
  !> set_ising_allup_potts: Set spins `1`.
  pure subroutine set_ising_allup_potts(this)
    class(potts), intent(inout) :: this
    this%spins_(:) = 1_int32
  end subroutine set_ising_allup_potts
  !> set_ising_allup_potts: Set spin to `1` in the first half of the region; in the other half, set it randomly.
  impure subroutine set_ising_mixed_phase_potts(this)
    class(potts), intent(inout) :: this
    integer(int64) :: half
    real(real64), allocatable :: r(:)
    !> ..... [16 ~ 20]
    !> ..... [11 ~ 15]
    !> ..... [ 6 ~ 10] 10 == (ny / 2) * nx
    !> ..... [ 1 ~  5]
    half = (this%ny_ / 2) * this%nx_
    this%spins_(1:half) = 1_int32
    allocate(r(half + 1:this%nall_))
    call random_number(r)
    this%spins_(half + 1:this%nall_) = floor(r * this%max_state_) + 1
    call this%update_norishiro()
  end subroutine set_ising_mixed_phase_potts
  !> set_ising_random_potts: Set spin to an integer between `1` and `this%max_state_` randomly.
  impure subroutine set_ising_random_potts(this)
    class(potts), intent(inout) :: this
    real(real64), allocatable :: r(:)
    allocate(r(1:this%nall_))
    call random_number(r)
    this%spins_(1:this%nall_) = floor(r * this%max_state_) + 1
    call this%update_norishiro()
  end subroutine set_ising_random_potts
  !> set_kbt_potts: Set parameter `beta` as `1 / kbt`.
  pure subroutine set_kbt_potts(this, kbt)
    class(potts), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_potts
  !> set_beta_potts: Set parameter `beta` and update `this%exparr_`.
  pure subroutine set_beta_potts(this, beta)
    class(potts), intent(inout) :: this
    real(real64), intent(in) :: beta
    integer(int32) :: center_before, center_after, i, j, k, l
    this%beta_ = beta
    do center_after = 1, this%max_state_
       do center_before = 1, this%max_state_
          do l = 1, this%max_state_
             do k = 1, this%max_state_
                do j = 1, this%max_state_
                   do i = 1, this%max_state_
                      associate(delta_e => - (count(center_after == [i, j, k, l]) - count(center_before == [i, j, k, l])))
                        if (delta_e <= 0.0_real64) then
                           this%ws_(i, j, k, l, center_before, center_after) = 1.0_real64
                        else
                           this%ws_(i, j, k, l, center_before, center_after) = exp(- this%beta_ * delta_e)
                        end if
                      end associate
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine set_beta_potts

  !> update_potts: Update the system by Metropolis method.
  impure subroutine update_potts(this)
    class(potts), intent(inout) :: this
    real(real64), allocatable :: r(:), next_state(:)
    integer(int64) :: i, j
    allocate(r(this%nall_), next_state(this%nall_))
    call random_number(r)
    do j = 1, 2
       do i = j, this%nall_, 2
          call this%update_onesite(i, r(i), floor(next_state(i) * this%max_state_) + 1)
       end do
       call this%update_norishiro()
    end do
  end subroutine update_potts
  !> update_onesite_potts: Update a spin of the system.
  pure subroutine update_onesite_potts(this, idx, r, next_state)
    class(potts), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r
    integer(int32), intent(in) :: next_state
    if (r < this%ws_( this%spins_(idx + 1), this%spins_(idx - 1), &
                    & this%spins_(idx + this%nx_), this%spins_(idx - this%nx_), &
                    & this%spins_(idx), next_state )) then
       this%spins_(idx) = next_state
    end if
  end subroutine update_onesite_potts
  !> update_norishiro_potts: Update norishiro.
  pure subroutine update_norishiro_potts(this)
    class(potts), intent(inout) :: this
    integer(int64) :: i
    do i = 1_int64, this%nx_
       this%spins_(this%norishiro_begin() + i - 1) = this%spins_(this%nall_ - this%nx_ + i)
       this%spins_(this%norishiro_end() - this%nx_ + i) = this%spins_(i)
    end do
  end subroutine update_norishiro_potts

  !> calc_total_energy_potts: Calculate the total energy.
  pure integer(int64) function calc_total_energy_potts(this) result(res)
    class(potts), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res - count(this%spins_(i) == [this%spins_(i + 1), this%spins_(i + this%nx_)])
    end do
  end function calc_total_energy_potts
  !> calc_total_magne_potts: Calculate the total magne.
  pure integer(int64) function calc_total_magne_potts(this) result(res)
    class(potts), intent(in) :: this
    integer(int64) :: i
    res = 0_int64
    do i = 1_int64, this%nall_
       res = res + merge(1, 0, this%spins_(i) == 1)
    end do
  end function calc_total_magne_potts

  !> nx_potts: Return size of `x` of the system.
  pure integer(int64) function nx_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%nx_
  end function nx_potts
  !> ny_potts: Return size of `y` of the system.
  pure integer(int64) function ny_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%ny_
  end function ny_potts
  !> nall_potts: Return size of the system.
  pure integer(int64) function nall_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%nall_
  end function nall_potts
  !> max_state_potts: Return maximum state of the system.
  pure integer(int64) function max_state_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%max_state_
  end function max_state_potts
  !> kbt_potts: Return temperature of the system.
  pure real(real64) function kbt_potts(this) result(res)
    class(potts), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_potts
  !> beta_potts: Return inverse temperature of the system.
  pure real(real64) function beta_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%beta_
  end function beta_potts
  !> norishiro_begin_potts: Return start index of `this%spins_(:)`.
  pure integer(int64) function norishiro_begin_potts(this) result(res)
    class(potts), intent(in) :: this
    res = 1 - this%nx_
  end function norishiro_begin_potts
  !> norishiro_end_potts: Return end index of `this%spins_(:)`.
  pure integer(int64) function norishiro_end_potts(this) result(res)
    class(potts), intent(in) :: this
    res = this%nall_ + this%nx_
  end function norishiro_end_potts
end module potts_m
