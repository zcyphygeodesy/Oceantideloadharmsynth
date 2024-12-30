## Fortran codes for spherical harmonic synthesis of ocean tidal load effects on all-element geodetic variations
https://www.zcyphygeodesy.com/en/h-nd-116.html
## [Algorithm purpose]
    Using the global ocean tidal load spherical harmonic coefficient model (cm), and given the longitude, latitude, orthometric height and forecast time of the calculation point, predict the sea surface tidal height (cm) and compute the ocean tidal load effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (10μE) or horizontal gravity gradient (NW, to the north and to the west, 10μE).
    Expand and improve the ocean tidal load effect algorithm in the IERS conventions (2010) to adapt to all-element geodetic variations in the whole Earth space.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAgsrbQuQYo4LKTTjCWDji6CQ.jpg)
## [Computation Output]
    tdh: the sea surface tidal height (cm). tdn(14): ocean tidal load effects on all-element geodetic variations.
    tdn(1:14) stores the ocean tidal load effects on 10 kinds of geodetic variations, which are the ocean tidal load effects on height anomaly tdn(1) (mm), ground gravity #tdn(2) (μGal), gravity disturbance tdn(3) (μGal), ground tilt #tdn(4:5) (SW, to the south and to the west, mas), vertical deflection tdn(6:7) (SW, to the south and to the west, mas), horizontal displacement #tdn(8:9) (EN, to the east and to the north, mm), ground radial displacement #tdn(10) (mm), ground normal or orthometric height #tdn(11) (mm), radial gravity gradient tdn(12 )(10μE) or horizontal gravity gradient tdn(13:14) (NW, to the north and to the west, 10μE).
    The calculation point can be on the ground, low altitude, satellite, ocean or underwater space. The geodetic variations abvove marked with # are valid only when the site is fixed with the solid Earth.
## [Geophysical models]
    (1) Global ocean tidal load spherical harmonic coefficient model FES2004S1.dat from IERS conventions 2010.
    (2) The Earth’s Load Love number file love_load_cm.dat from a Regional EIAstic Rebound calculator (REAR1.0, 2015).
## [Main program for test entrance]
    Otdloalharmsynth.f90
    The record of the test output file reslt.txt: the difference between the MJD day and starting MJD0, tdh, tdn(1:14)
    The test program reads two geophysical models into memory in advance, The calculation time for one calculation point is about 25 ms when the maximum trucated degree set as 120.
## (1) Algorithm module for global ocean tidal height forecast by spherical harmonic synthesis
    otidehsynth(mjd,BLH,tdh,fes,nn,maxn,GRS)
    Input parameters: mjd – the epoch time in MJD.
    Input parameters: BLH(3) - latitude, longitude (decimal degrees) and orthometric height (m) at the calculation point.
    Input parameters: fes(nn,7) - the tidal load spherical harmonic coefficients [doodson,n,m,C+(cm),eps+,C-(cm),eps-].
    Input parameters: maxn - maximum trucated degree.
    Input parameters: GRS(6) - gm,ae,j2,omega,1/f, default value.
    Return parameters: tdh - the sea surface tidal height (cm).
## (2) Algorithm module for spherical harmonic synthesis of the load effects on  all-element geodetic variations
    LTideFlupnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr)
    Input parameters: rln(3) - the spherical coordinates of calculation point in IERS.
    Input parameters: flv(maxn,3) - the load Love numbers.
    Input parameters: pnm,dpt1,dpt2,gr - the normalized associative Legendre functions pnm and their derivatives dpt1,dpt2 and normal gravity at the calculation point.
    Return parameters: tdn(1:14) - the ocean tidal load effects on all-element geodetic variations.
## (3) Algorithm module for the direct influence of ocean tidal load to geopotential coefficient variations
    OLoadDFlu(mjd,cnm,snm,maxn,fes,nn)
    Return parameters: cnm,snm - the direct influence of ocean tidal load to geopotential coefficient variations.
## (4) Algorithm module for normalized associative Legendre functions and their derivatives
    BelPnmdt(pnm,dpt1,dpt2,maxn,t)
    Improved Belikov recursion algorithm for pnm and non-singular recursive algorithm for derivative of pnm.
## (5) Algorithm module for the nodal corrections of the tidal constituent
    CalcTidefu(mjd,doodson,df,du)
    Input parameters: doodson(6) - Doodson constants of the tidal constituent.
    Return parameters: df,du - the nodal factor and nodal correction angle (degree) of the tidal constituent.
## (6) Algorithm module for the phase bias of the tidal constituent
    BiasTide(doodson,bias)
    Return parameters: bias - the phase bias (degree) of the tidal constituent.
## (7) Algorithm module for the 5 basic astronomical mean longitudes
    ASTRO5(TIME,SHPNP) ! s, h, p, N, p'
## (8) Calculation module for the modified Julian Date
    CAL2JD(IY,IM,ID,mjd,j)  !From IAU SOFA library
## (9) Algorithm module for transforming ellipsoid geodetic coordinates into spherical coordinates
    BLH_RLAT(GRS,BLH,RLAT)
## (10) Calculation module for normal geopotential coefficients
    normdjn(GRS,djn) !GRS(6) - gm,ae,j2,omega,1/f, default value
## (11) Calculation module for Legendre functions and their derivatives to ψ
    Legendre(maxn,m,rlat,RLEG,DLEG); PlmBar_d(p,dp,lmax,rlat)
## (12) Algorithm module for transforming the long integer time (date) agreed by ETideLoad into year, month, day, hour, minute and second.
    tmcnt(tm,iyr,imo,idy,ihr,imn,sec)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] ETideLoad4.5 User Reference https://www.zcyphygeodesy.com/en/
    8.4.1 Construction of tidal load spherical harmonic coefficient model
    8.2.2 The normalized spherical harmonic series expansion for surface load deformation field
    8.2.3 The normalized associated Legendre functions and thier derivatives
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file, geophysical models and all input and output data.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAgsrbQuQYopIaw9QQwlg44ugk.jpg)
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAgs7bQuQYonJOz9QMwlg44ugk.jpg)
