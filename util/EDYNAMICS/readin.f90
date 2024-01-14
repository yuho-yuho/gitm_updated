
module readhlpoten_module 

  ! Global
  use params_module, only: nmlon, nmlat_T1
      
  implicit none
 
  ! necessary parameters for calc_FAC                                     
  real, dimension(nmlon,nmlat_T1) :: poten_hl       

end module readhlpoten_module
      
