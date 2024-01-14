program test_efvm

  use efvm

  implicit none    

  real :: er, disp, scale, scale_sh
                                                                               
  write(*,*) "Read IMF By and Bz"                                              
  read(*,*) IO_IMFBy, IO_IMFBz                                                 
                                                                               
  IO_SWVX=450.                                              
  IO_SWN=4.                        

  er=1.2
  disp=0.
  scale=1.4
  scale_sh=1.4
                                                                               
  call initialize_efvm                                 
  call efvm_main(er,disp,scale,scale_sh)

  if (efvm_debug) then
     write(*,*) "------ Check outputs ------"
     write(*,*) maxval(efvm_dEd1), maxval(efvm_dEd2)
     write(*,*) maxval(efvm_dEd1_sh), maxval(efvm_dEd2_sh)
  end if
end program test_efvm
