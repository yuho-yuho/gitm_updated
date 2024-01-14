! This Module serves as a new electrodynamics module for GITM
! Created: Qingyu Zhu, 07/26/17, qingyu.zhu@mavs.uta.edu
!
! Using part of the new NCAR 3D dynamo module (Maute and Richmond, SSR, 2016)
!
! Currently, the calculation is not parallized 
!
! The module is driven by: 
!        1. Zonal and horizontal winds
!        2. Pedersen and Hall conductivities 
!        3. High latitude Potential => Field-aligned current (RHS of eqs)
! These parameters are collected from each procedure by MPI 
! 
! The parameters output from the module are:
!         mlat, mlon, potential, fac
! These parameters are then spread to each procedure
!
! Add in this GITM on 03/02/2020, qingyu
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ModGedy

  use fieldline_p_module,only: fieldline_p_dim
  use fieldline_s_module,only: fieldline_s_dim
  use fieldline_r_module,only: fieldline_r_dim
  use qd_module,only: gen_qd_grid

  use coef_module,only: calc_coef,add_coef_ns,const_rhs, &
       calc_FAC, calc_coef1, correct_potential, close_winddynamo

  use area_factors_module, only: calc_mf
  
  use delB_module, only: gen_geo_grid

  ! Need to specify what output parameters later
  use edynamicsout_module

  ! FAC driven modules 
  ! qingyu, 11/24/2020
  ! Module(s) from the EDYNAMICDIR
  !use ModFAC, only: const_init_coef_ns2, get_asym_pot, remove_jwd
  use ModFAC, only: get_asym_pot

  ! Module(s) from current directory
  use ModAMPERE, only: read_ampere, get_currenttime_fac
  use ModOIM, only: read_oimfac, get_currenttime_oimfac
  use ModSDFM, only: read_sdfm, get_currenttime_sdfm

  use ModInputs, only: UseAMPERE, UseASYFAC, UseOIM,  &
       UseWindDynamo, RemoveJwd, UseSDFM

  implicit none

  ! Global variables   
  !!!!!!! GITM variables
  ! Wind and conductivities
  real, allocatable :: vn_gitm(:,:,:), un_gitm(:,:,:), &
       sigP_gitm(:,:,:), sigH_gitm(:,:,:)

  ! 2D Potential in geomagnetic coordinate (for high latitude potential)
  ! Also output 2D ed1 and ed2, Qingyu, 01/24/2018
  integer, parameter :: nmlon_gitm=144, nmlat_gitm=72
  real, allocatable :: mlon1_gitm(:), mlt1_gitm(:), mlat1_gitm(:)
  real, allocatable :: poten2d_gitm(:,:)
  real, allocatable :: ed1_2d_gitm(:,:), ed2_2d_gitm(:,:)

  ! Geographic coordinates (mapping)
  real, allocatable :: lon_gitm(:), lat_gitm(:), hgt_gitm(:)

  ! Magnetic coordinates (interpolation)
  real, allocatable :: mlon_gitm(:,:,:,:), mlat_gitm(:,:,:,:)
  integer :: nlon_gitm, nlat_gitm, nhgt_gitm

  ! Build the structure to store the info for re-arrangement
  type task
     real :: lon0 ! starting lon (left)
     real :: lat0 ! starting lat (bottom)
  end type task

  type(task), allocatable :: tasks(:)

  contains

    !************************************************************************
    
    subroutine init_gedy(year,month,day)

      ! This routine initializes the 3D dynamo model
      ! Created: Qingyu Zhu, 07/26/2017

      implicit none

      ! Local
      integer, intent(in) :: year, month, day
      integer :: doy
      logical, parameter :: debug =.true.
      real :: date_run

      ! Set up magnetic latitude and longitude grid                          
      call gen_highres_grid
      if(debug) write(6,*) '=> Done gen_highres_grid!'

      ! Set up p-fieldline grid                                       
      call fieldline_p_dim
      if(debug) write(6,*) '=> Done fieldline_p_dim!'

      ! Set up s-fieldline grid                                            
      call fieldline_s_dim
      if(debug) write(6,*) '=> Done fieldline_s_dim!'

      ! Set up r-fieldline grid                                  
      call fieldline_r_dim
      if(debug) write(6,*) '=> Done fieldline_r_dim!'

      ! Set up qd grid
      call gen_qd_grid
      call gen_geo_grid
      if(debug) write(6,*) '=> Done gen_qd_grid!'

      ! get APEX
      call calc_doy(year,month,day,doy)
      date_run = year + doy/365.
      write(*,*) date_run
      call apxparm(date_run)
      if(debug) write(6,*) '=> Done apxparm!'

      ! Initialize AMPERE module
      ! qingyu, 11/24/2020
      if (UseAMPERE) then
         call read_ampere
         write(*,*) "=> Done initializing AMPERE"
      else if (UseOIM) then
         call read_oimfac
         write(*,*) "=> Done initializing OIM FAC"
      else if (UseSDFM) then
         call read_sdfm
         write(*,*) "=> Done initializing SDFM FAC"       
      end if

    end subroutine init_gedy
    
    !------------------------------------------------------------------------

    subroutine calc_gedy

      ! Routine to calculate the potential/fac 
      ! Need winds, conductivities and high latitude potential
      ! Created: Qingyu Zhu, 07/27/2017

      ! Local
      implicit none
      logical, parameter :: debug = .true.
      
      ! calculate M- and N-coefficients                                     
      !(S1,S2,R points, but for QD is done in get_apex)             
      call calc_mf
      call calc_mn_s1s2
      if(debug) write(6,*) '=> Done calc_MN_S1S2!'
      
      ! calculate JeD-coefficients (right hand side)            
      call calc_je_s1s2
      if(debug) write(6,*) '=> Done calc_je_s1s2!'
      
      ! calculate S (right hand side)                       
      call calc_S
      if(debug) write(6,*) '=> Done calc_S!'
      
      ! calculate coefficients (left hand side)           
      call calc_coef
      if(debug) write(6,*) '=> Done calc_coef!'

      ! Turn off the wind dynamo
      if (.not. UseWindDynamo) call close_winddynamo
    
      ! calculate high latitude forcing (left hand side)
      ! qingyu, 11/24/2020
      ! add the capability of using AMPERE FAC
      if (UseAMPERE) then
         call get_currenttime_fac
         call get_ampere_fac
      else if (UseOIM) then
         call get_currenttime_oimfac                         
         call get_oim_fac
      else if (UseSDFM) then
         call get_currenttime_sdfm                         
         call get_sdfm_fac
      end if
      
      if(debug) write(6,*) '=> Done Get FAC'

      ! qingyu, 11/25/2020
      ! For asymmetric FAC
      if ((UseAMPERE .or. UseOIM .or. UseSDFM) .and. UseASYFAC) then
         call get_asym_pot
         if(debug) write(6,*) '=> Done calc asymmetric FAC potential!'
      end if 

      write(6,*) "==> Done determine high-latitude electric potential"

      ! Determine the transition zone and modify the coefficients
      call calc_coef1 
      if(debug) write(6,*) '=> Done modify coef_ns2' 

      ! add NH+SH and sum up in height: lhs & rhs                             
      call add_coef_ns
      if(debug) write(6,*) '=> Done add coef_ns2 in different hemispheres!'
      
      ! construct LHS & RHS matrix and solve                                  
      call const_rhs
      if(debug) write(6,*) '=> Done calculate symmetric potential!'

      ! correct high-latitude electric potential since they do not 
      ! need to be identical 
      call correct_potential
      if(debug) write(6,*) '=> Done correct potential!'

      ! calculate electric fields (This is not time consuming)
      call calc_3d_Efield
      if(debug) write(6,*) '=> Done calc Efield from 3D model!'

      ! output the the results from fieldlines to matrixs                
      call edynamics_output
      if(debug) write(6,*) '=> Done output from 3D model!'
      
    end subroutine calc_gedy

end module ModGedy

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! A routine to calculate the DOY 

subroutine calc_doy(year,month,day,doy)

  implicit none

  integer, intent(in) :: year, month, day
  integer, intent(out) :: doy
  integer :: imon
  integer, parameter :: mo(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  logical :: leap=.false.

  doy=0

  do iMon=1,month-1
     doy=doy+mo(iMon)
  enddo

  doy=doy+day

  if (mod(year,100) .ne. 0 .and. mod(year,4) .eq. 0) leap = .true.
  if (mod(year,400) .eq. 0) leap = .true.

  if(leap) doy = doy +1

end subroutine calc_doy
