program test_epm

  use EPM_Interface                                                           
  use EPM_Initialization                                                      
  use epm                                                                     
                                                                               
  implicit none                                                                
                                                                               
  write(*,*) "Read IMF By and Bz"
  read(*,*) IO_IMFBy, IO_IMFBz
  
  IO_SWVX=450.                                               
  IO_SWN=4.                                                                    
                                                                               
  call initialize_epm                                                         
  call epm_main 

end program test_epm
