
program GITM

  use ModInputs
  use ModTime
  use ModGITM

  implicit none

  integer :: iBlock

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  call init_mpi
  call start_timing("GITM")
  call delete_stop

  call init_planet
  call set_defaults

  call read_inputs(cInputFile)
  call set_inputs

  call initialize_gitm(CurrentTime)

  call write_output

  call report("Starting Main Time Loop",0)

  ! ------------------------------------------------------------------------
  ! Run for a few iterations
  ! ------------------------------------------------------------------------

  do while (CurrentTime < EndTime)

     call calc_pressure

     !!! We may have to split cMax and Dt calculation!!!
     Dt = 1.e32

     call calc_timestep_vertical
     if (.not. Is1D) call calc_timestep_horizontal

     call advance

     ! Advance the Gedy, qingyu, 03/02/2020
     call advance_gedy

     ! Post-processes 
     ! Calculate the TEC, qingyu, 10/15/2020
     ! Calculate the height-integrated Joule heating, 11/26/2020
     do iBlock=1,nBlocks 
        call calc_tec(iBlock)
        call calc_integrated_jh(iBlock)
        call calc_on2ratio(iBlock)
     end do

     if (.not.IsFramework) call check_stop

     iStep = iStep + 1

     call write_output

  enddo

  ! ------------------------------------------------------------------------
  ! Finish run
  ! ------------------------------------------------------------------------

  call finalize_gitm

end program GITM

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  call stop_gitm(StringError)

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  implicit none
  integer, intent(in) :: iUnit

  return

end subroutine CON_io_unit_new

!---------------------------------------------------------------------------

