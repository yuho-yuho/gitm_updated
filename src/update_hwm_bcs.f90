! -----------------------------------------------------------------------------
! Update the lower boundary winds, introduced by Cissi Lin since 2016
! qingyu zhu, 03/02/2020
! -----------------------------------------------------------------------------

subroutine update_hwm_bcs

    use ModGitm
    use ModTime
    use ModInputs
    use ModKind, ONLY: Real8_
    use ModConstants, only: pi

    implicit none

    integer  :: iBlock, iSpecies, iLat, iLon, iAlt, iyd

    real     :: geo_lat, geo_lon, geo_alt, geo_lst
    real*4   :: hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst
    real*4   :: hwm_f107a, hwm_f107, hwm_ap(2), qw(2)


    iyd = iTimeArray(1)*1000 + iJulianDay


    do iBlock = 1, nBlocks
        do iLat = 1, nLats
            do iLon = 1, nLons
                do iAlt = 1, 2

                    geo_lat = Latitude(iLat,iBlock)*180.0/pi
                    geo_lon = Longitude(iLon,iBlock)*180.0/pi

                    geo_alt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
                    geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)

                    hwm_utime = utime
                    hwm_alt = geo_alt
                    hwm_lat = geo_lat
                    hwm_lon = geo_lon
                    hwm_lst = geo_lst
                    hwm_f107a = f107a
                    hwm_f107 = f107
                    hwm_ap(1) = -1.0
                    hwm_ap(2) = -1.0
              
                    call HWM14(iyd,hwm_utime,hwm_alt,hwm_lat,hwm_lon,hwm_lst,&
                         hwm_f107a,hwm_f107,hwm_ap,qw)

                    ! qw is north&east
                    Velocity(iLon,iLat,iAlt,iEast_,iBlock) = qw(2)
                    Velocity(iLon,iLat,iAlt,iNorth_,iBlock) = qw(1)

                    ! set vertical wind to 0 at LB
                    VerticalVelocity(iLon, iLat, iAlt, :, iBlock) = 0.

                enddo
            enddo
        enddo
    enddo 


end subroutine update_hwm_bcs
