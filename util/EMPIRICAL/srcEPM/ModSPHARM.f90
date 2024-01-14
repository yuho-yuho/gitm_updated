! This routine serves as the module for the SPHARM fitting 
!
! Created: Qingyu Zhu, 11/03/2019
!
! ----------------------------------------------------------------------------

module modspharm

  implicit none

  contains

    ! ------------------------------------------------------------------------
    ! Caculate the Associated Legengre Functions (ALF) by using the 
    ! recursion relationship
    ! ALFs: P(l,m): l>0, m>0

    subroutine calc_alf(x,LMAX,alfs)

      ! LMAX: Fitting order
      ! x: cos colat
      ! a matrix with the size of (LMAX+1) x (LMAX+1)

      implicit none 

      real, intent(in) :: x
      integer, intent(in) :: LMAX
      real, intent(out) :: alfs(LMAX+1,LMAX+1)

      integer :: l, il, m, im,  dfactorial
      real :: val1, val2

      alfs=0.

      alfs(1,1) = 1.

      do il=2,LMAX+1

         l=il-1

         do m=0,l-1

            im=m+1
            
            ! P(l-1,m)
            if (l-1<m) then
               val1=0.
            else
               val1=alfs(il-1,im)
            end if

            ! P(l-2,m)
            if (l-2<m) then
               val2=0.
            else
               val2=alfs(il-2,im)
            end if
               
            alfs(il,im) = (x*(2*l-1)*val1 - (l+m-1)*val2) / (l-m)

         end do
 
         ! Calculate P(l,l)
         alfs(il,il)=((-1.)**(l))*dfactorial(2*l-1)*(1-x**2)**(l*0.5)

      end do

    end subroutine calc_alf

    ! --------------------------------------------------------------------
    ! Reconstruct the cos(m*phi)P(l,m) and sin (m*phi)P(l,m)
    ! for a certain location
    ! Works only for the Coeffs for our model due to the order

    subroutine reconst_spharm(theta,phi,LMAX,spharms) 
    
      ! theta (colat) --> x = cos(theta), 0<= theta <= pi
      ! phi, 0<= phi <= 2*pi

      implicit none

      real, intent(in) :: theta, phi
      integer, intent(in) :: LMAX
      real, intent(out) :: spharms((LMAX+1)**2)

      integer :: l, il, m, im, loc
      real :: x, alf, mp, cosmp, sinmp
      real :: alfs(LMAX+1,LMAX+1)

      ! First calculate the Associated Legendre Functions at x
      x=cos(theta)
      call calc_alf(x,LMAX,alfs)
      
      !Construct the cos(m*phi)P(l,m) and sin(m*phi)P(l,m) for all l and m
      spharms=0.

      ! First reconstruct the case m=0 for all l
      do il = 1,LMAX+1
         spharms(il)=alfs(il,1)
      end do

      ! Then reconstruct for each l and m<=l
      loc = LMAX + 2
      
      do l=1,LMAX

         il=l+1

         do m=1,l

            ! I think that it might be better to keep this
            
            !if (m>3) cycle

            im=m+1
            alf=alfs(il,im)
            mp=m*phi

            cosmp=cos(mp)
            sinmp=sin(mp)

            spharms(loc)=cosmp*alf
            loc=loc+1

            spharms(loc)=sinmp*alf
            loc=loc+1
         end do
         
      end do

      !write (*,*) spharms

    end subroutine reconst_spharm

end module modspharm

! ------------------------------------------------------------------------
! Double Factorial

recursive function dfactorial(n) result(fact)                 
  
  implicit none                                         
  
  integer :: fact                                                     
  integer, intent(in) :: n                                     
  
  if (n<=1) then                                     
     fact = 1                                                
  else                                                              
     fact = n * dfactorial(n-2)                            
  end if
  
end function dfactorial

