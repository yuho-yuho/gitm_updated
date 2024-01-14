! New routine used to update the N2A by Cissi Lin
! Added in this version by qingyu, 03/03/2020

subroutine update_n2a(iBlock)

    use ModEUV, only: Sza
    use ModGITM, only: nLats, nLons, nAlts, nBlocks, NDensityS, Altitude_GB, Latitude, Longitude
    use ModInputs
    use ModNOChemistry
    use ModPlanet, only: iN2_A_
    use ModTime

    implicit none

    integer, intent(in)      :: iBlock

    character(len=60)        :: N2A_File
    real, dimension(24)      :: data_line
    integer                  :: iLon, iLat, iAlt, geo_lst, ii, idx, idx_next
    integer, dimension(1:7)  :: itime
    real                     :: minsec, alt, geo_lon
    real                     :: N2A_interp, N2A_interp_HrNow, N2A_interp_HrNext

    !itime(1) = nyear + 1965
    !itime(2) = nmonth
    !itime(3) = nday + 1
    !itime(4) = nhour
    !itime(5) = nmin
    !itime(6) = nsec
    !itime(7) = (timeleft - nsec) * 1000
    call time_real_to_int(CurrentTime, itime)

    call report("update_n2a",1)
    call start_timing("update_n2a")

    if (iDebugLevel > 0) write(*,*) "==> # of times in N2(A) : ", nN2A_times

    ! read in N2A files
    if (.not. IsN2Aread) then

        IsN2Aread = .true.

        ! prepare altitudinal grids (in km) for N2A profiles
        N2A_ZZ = (/ (ii, ii = 1, 106) /)
        N2A_ZZ = N2A_ZZ * 2.0 + 38.


        ! File name of N2A Profiles
        N2A_File ='UA/DataIn/N2A/n2a.csv'

        ! open and read the data file
        if (iDebugLevel > 0) write(*,*) "==> Reading File : ", N2A_File

        open(unit=128, file=N2A_File, status='old', action='read')

        do ii = 1, 106

            read(128, *) data_line
            N2A(ii,:) = data_line * 1e6

        enddo

        close(unit=128)

    endif



    ! calculate N2A profile at the location
    do iLon = -1, nLons+2
        do iLat = -1, nLats+2
            do iAlt = -1, nAlts+2

                ! calculate local sunlit time
                geo_lon = Longitude(iLon,iBlock) * 180.0 / pi
                geo_lst = mod(utime/3600.0 + geo_lon / 15.0, 24.0)


                if (cos(sza(iLon, iLat, iBlock)) > 0.0) then

                    ! interpolate N2(A) profiles
                    minsec = itime(5) * 60. + itime(6) 
                    alt = Altitude_GB(iLon, iLat, iAlt, iBlock) / 1000.0
                    idx = minloc(abs(N2A_ZZ - alt), 1)
                    if (N2A_ZZ(idx) >= alt) idx_next = idx - 1
                    if (N2A_ZZ(idx) < alt) idx_next = idx + 1
                    ! #1 - interpolate in altitude at the same hours
                    N2A_interp_HrNow = (N2A(idx, geo_lst) * abs(N2A_ZZ(idx_next) - alt) &
                                   + N2A(idx_next, geo_lst) * abs(N2A_ZZ(idx) - alt)) / 2.0
                    N2A_interp_HrNext = (N2A(idx, geo_lst+1) * abs(N2A_ZZ(idx_next) - alt) &
                                    + N2A(idx_next, geo_lst+1) * abs(N2A_ZZ(idx) - alt)) / 2.0
                    ! #2 - interpolate in time
                    N2A_interp = (N2A_interp_HrNow * (3600.0 - minsec) + N2A_interp_HrNext * minsec) / 3600.0

                    NDensityS(iLon, iLat, iAlt, iN2_A_, iBlock) = N2AScaling * cos(sza(iLon, iLat, iBlock)) * N2A_interp &
                      * (f107 + f107a) / 2.0 / 180.0  ! this line scales the N2(A) density which is generated at P=180

                    if ((NDensityS(iLon, iLat, iAlt, iN2_A_, iBlock) < 0.0)) then
                    write(*,*) "Negative N2(A) Density!!!"
                    write(*,*) "idx, idx_next = ", idx, idx_next
                    write(*,*) "N2A_HrNow = ", N2A(idx, geo_lst), N2A(idx_next, geo_lst)
                    write(*,*) "N2A_HrNext = ", N2A(idx, geo_lst+1), N2A(idx_next, geo_lst+1)
                    write(*,*) "abs(height difference) = ", abs(N2A_ZZ(idx_next) - alt), abs(N2A_ZZ(idx) - alt)
                    write(*,*) "N2A_interp = ", N2A_interp_HrNow, N2A_interp_HrNext, N2A_interp
                    write(*,*) "cos = ", cos(sza(iLon, iLat, iBlock))
                    write(*,*) "geo_lst = ", geo_lst
                    write(*,*) "N2(A) = ", NDensityS(iLon, iLat, iAlt, iN2_A_, iBlock)
                    endif

                else
                    NDensityS(iLon, iLat, iAlt, iN2_A_, iBlock) = 0.0
                endif

            enddo
        enddo
    enddo

    if (iDebugLevel > 0) write(*,*) "Finished updating N2(A)!"


    call end_timing("update_n2a")


end subroutine update_n2a
