! NO chemistry module developed by Cissi Lin
! Added in this GITM by qingyu, 03/03/2020

module ModNOChemistry

    implicit none

    logical                 :: IsN2Aread = .false.
    logical                 :: UseNOHigherVibrationalStates
    logical                 :: UseTCIParameterizationS, UseTCIParameterizationGM
    real, dimension(106,24) :: N2A
    real, dimension(106)    :: N2A_ZZ
    real                    :: N2AScaling, ScalingIonPrecipIonRate
    real                    :: ScalingEuvAbsRate, ScalingEuvIonRate
    real                    :: VSP_f = 15.0
    integer                 :: nN2A_times
    real, dimension(10)     :: Adv1, Adv2, VSP_k
    real, dimension(11)     :: VSP_g2d, VSP_g4s


    data Adv1 /12.396, 17.901, 22.974, 27.637, 31.911, &
               35.816, 39.374, 42.604, 45.523, 48.150 /

    data Adv2 /0.0, 0.786, 1.534, 2.489, 3.630, &
               4.932, 6.372, 7.927, 9.574, 11.292/

    ! k in unit of m^3/s
    data VSP_k /42.0e-18, 12.390e-18, 8.7567e-18, 6.938e-18, 5.531e-18, &
                5.101e-18, 4.302e-18, 4.011e-18, 3.578e-18, 3.449e-18/

    data VSP_g2d /0.0552162, 0.0520401, 0.0508185, 0.0408014, &
                  0.0522844, 0.0752504, 0.111654,  0.135353,  &
                  0.144149,  0.121915,  0.160518/
    data VSP_g4s /0.106307,  0.116395,  0.117199, 0.147256,  &
                  0.168128,  0.127784,  0.101884, 0.0682910, &
                  0.0292211, 0.0122811, 0.00525504/

end module ModNOChemistry
