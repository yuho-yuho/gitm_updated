! Solver (qingyu, 11/26/2020)
! Adapted from subroutine const_rhs
! -------------------------------------------------------------------------
subroutine pot_solver(coef_in,pot_out)

  use params_module, only: nmlon, nmlat_h, nlonlat 
  use fieldline_p_module, only: fline_p
  use coef_module, only: lhs_ns, rhs_ns

  implicit none
  
  real, intent(in) :: coef_in(nmlon,nmlat_h,10)
  real, intent(out) :: pot_out(nmlon,nmlat_h)
      
  integer :: i,j,im,ip,ijm,it,ij,jt,istat,isn

  integer:: info,ifail
  integer,parameter :: nrhmax = 1
  integer ::  ipiv(nlonlat)
  character, parameter :: trans='N'
  
  real(8) :: colsum,anorm ,rcond,work(4*nlonlat) 
  real(8) :: dlange 
  integer :: info_c,iwork(nlonlat)
  character, parameter :: norm='1' 
  
  if (.not. allocated(lhs_ns)) then
     allocate(lhs_ns(nlonlat,nlonlat),STAT=istat)
     allocate(rhs_ns(nlonlat),STAT=istat)
  end if
  
  lhs_ns=0.
  rhs_ns=0.
  ij=0
  isn=1
  pot_out=0.
  
  do i=1,nmlon
     
     if (i==1) then
        im=nmlon
     else
        im=i-1
     end if
     
     if (i==nmlon) then
        ip=1
     else
        ip=i+1
     end if
     
     j=1
     ij = (i-1)*nmlat_h+j
     rhs_ns(ij) = coef_in(i,j,10)
     
     it = (ip-1)*nmlat_h+j
     lhs_ns(ij,it)  = coef_in(i,j,1)
     
     it = (ip-1)*nmlat_h+j+1 
     lhs_ns(ij,it)  = coef_in(i,j,2)
     
     it = (i-1)*nmlat_h+j+1
     lhs_ns(ij,it)  = coef_in(i,j,3)
     
     it = (im-1)*nmlat_h+j+1
     lhs_ns(ij,it)  = coef_in(i,j,4)
     
     it = (im-1)*nmlat_h+j
     lhs_ns(ij,it)  = coef_in(i,j,5)
     
     it = (i-1)*nmlat_h+j 
     lhs_ns(ij,it)  = coef_in(i,j,9)
     
     do j=2,nmlat_h-1
        
        ij = (i-1)*nmlat_h+j                                               
        rhs_ns(ij) = coef_in(i,j,10)                            
        
        it = (ip-1)*nmlat_h+j                               
        lhs_ns(ij,it)  = coef_in(i,j,1)             
        
        it = (ip-1)*nmlat_h+j+1                                            
        lhs_ns(ij,it)  = coef_in(i,j,2)                
        
        it = (i-1)*nmlat_h+j+1                                             
        lhs_ns(ij,it)  = coef_in(i,j,3)
            
        it = (im-1)*nmlat_h+j+1             
        lhs_ns(ij,it)  = coef_in(i,j,4)     
        
        it = (im-1)*nmlat_h+j                                              
        lhs_ns(ij,it)  = coef_in(i,j,5)                           
        
        it = (im-1)*nmlat_h+j-1                              
        lhs_ns(ij,it)  = coef_in(i,j,6)                          
        
        it = (i-1)*nmlat_h+j-1                                
        lhs_ns(ij,it)  = coef_in(i,j,7)                   
        
        it = (ip-1)*nmlat_h+j-1                         
        lhs_ns(ij,it)  = coef_in(i,j,8)                      
        
        it = (i-1)*nmlat_h+j                                            
        lhs_ns(ij,it)  = coef_in(i,j,9)
     end do
     
     j=nmlat_h                                                       
     ij = (i-1)*nmlat_h+j                                      
     rhs_ns(ij) = coef_in(i,j,10)                                
     
     it = (ip-1)*nmlat_h+j                                       
     lhs_ns(ij,it)  = coef_in(i,j,1)                                
     
     it = (im-1)*nmlat_h+j                                       
     lhs_ns(ij,it)  = coef_in(i,j,5)                              

     it = (im-1)*nmlat_h+j-1                                      
     lhs_ns(ij,it)  = coef_in(i,j,6)                        
     
     it = (i-1)*nmlat_h+j-1                                      
     lhs_ns(ij,it)  = coef_in(i,j,7)                       
     
     it = (ip-1)*nmlat_h+j-1                                  
     lhs_ns(ij,it)  = coef_in(i,j,8)                          
     
     it = (i-1)*nmlat_h+j                    
     lhs_ns(ij,it)  = coef_in(i,j,9)
     
  end do

  ! Solve the LHS X = RHS
  anorm = 0.d0
  
  do i=1,nlonlat
     colsum=0.d0
     
     ! Compute norm-1 of A
     do j=1,nlonlat
        colsum = colsum + abs(lhs_ns(j,i))
     end do
     
     anorm = max(anorm,colsum)
     
  end do
  write(6,*) 'anorm ',anorm
  
  write (6,*) 'factorize A'
  call DGETRF(nlonlat,nlonlat,lhs_ns,nlonlat,ipiv,info)
  write(6,*) 'info ' , info
  
  if (info .ne. 0) then
     write(*,*) "Singular matrix"
     write(*,*) "Error in LHS and RHS constructions, check !!!"
  else
     CALL DGECON(norm,nlonlat,lhs_ns,nlonlat,anorm,rcond,&
          work,iwork,info_c)
     write(6,*) 'info_c ' , info_c 
     write(6,*) 'condit=', 1/rcond
     
     write (6,*) 'solve'
     call DGETRS(trans,nlonlat,nrhmax,lhs_ns,nlonlat,ipiv,&
          rhs_ns,nlonlat,info)
     write (6,*) 'after solve',info 

     it=0
     do i=1,nmlon
        do j=1,nmlat_h
           
           it=it+1
           pot_out(i,j)=rhs_ns(it)

        end do
     end do
     
  end if
  
end subroutine pot_solver

