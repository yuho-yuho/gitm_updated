! To advance the New electrodynamics 
! update every few minutes
! Created: Qingyu Zhu, qingyu.zhu@mavs.uta.edu
! Add in this GITM on 03/02/2020, qingyu  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine advance_gedy

  ! Global
  use ModMpi
  use ModGITM
  use ModInputs
  use ModGedy
  use ModTime

  implicit none

  ! Local
  integer :: i, iBlock, iLon, iLat, iAlt, iError
  integer, parameter :: root = 0
  integer :: istep_tie = 1, istep_end = 1

  ! ---------------------------------------------------------------------------

  if (useGedy) then

     call report("advance_gedy",1)
     call start_timing("advance_gedy")

     if (floor((tSimulation-dt)/DtGedy) /= &
          floor((tsimulation)/DtGedy) .or.isFirstGedy) then

        !write(*,*)
        !write(*,*) '****** Updating Gedy ******'
        !write(*,*)

        !!!!!!!! First, Gathering different fields

        ! MPI gather (un, vn, sigP, sigH)
        call gedy_gatherfield

        ! MPI gather mlon, mlat
        call gedy_sendmagcoordinate
        
        if (iProc==0) then

           ! Generate 2D potential distribution in the root processor
           call gedy_get2dpoten

           ! Map fields to S1 and S2 points 
           call gedy_mapfield_s1s2

           ! Map fields to P points
           call gedy_mapfield_p

           ! Caculate the Electrodynamics in Gedy
           call calc_gedy

        endif

        ! Spread field from root and get 3D values                            
        call gedy_spreadfield_new
        call gedy_get3dvalues_new
        if (iProc==0) write(*,*) "=>Done spread fileds"

        isFirstGedy = .false.

     endif

     call end_timing("advance_gedy")

  endif

end subroutine advance_gedy
