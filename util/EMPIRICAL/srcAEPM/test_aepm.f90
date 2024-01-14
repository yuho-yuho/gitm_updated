program test_aepm

  use AEPM_Interface
  use AEPM_Initialization
  use aepm

  implicit none 

  write (*,*) "Hello World"

  read(*,*) IO_IMFBy, IO_IMFBz 
  
  IO_SWVX=450.
  IO_SWN=4.

  call initialize_aepm
  call aepm_main



end program test_aepm
