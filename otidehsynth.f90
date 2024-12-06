      subroutine otidehsynth(mjd,BLH,tdh,fes,nn,maxn,GRS)
	integer::nn,i,j,n,m,k,maxn
	integer*4::bs,td
	real*8::tdt,tdh,fes((maxn+2)**2*40,7),BLH(3),rln(3),rr,rlat,rlon,t,pi,RAD,bias
	real*8::sigma,th1,th2,astr(6),tt,dod(6),nm,cp,cn,sp,sn,ksp,ksn
	real*8::pnm(900000),pnm1(900000),mjd,df,du,GRS(6)
!---------------------------------------------------------------------
      call BLH_RLAT(GRS,BLH,rln)
	pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      rr=rln(1);rlat=rln(2);rlon=rln(3)*RAD
      call PlmBar_d(pnm,pnm1,maxn+1,rlat)
 	t=(mjd-51544.5d0)/36525.d0 !J2000.0 世纪
	tt=(mjd-51545.d0)-floor(mjd-51545.d0)!补回0.5天
      astr(2)=218.31664563d0+481267.88119575d0*t-0.001466388889d0*t**2    !s
     >	-0.000000074112d0*t**3-0.000000153389d0*t**4
	astr(3)=280.4664501606d0+36000.769748805556d0*t+0.000303222222d0*t**2 !h
     >	-0.000001905501d0*t**3-0.000000065361d0*t**4
	astr(4)=83.35324312d0+4069.01363525d0*t-0.010321722222d0*t**2     !p
     >	-0.000014417168d0*t**3+0.0000000526333d0*t**4
	astr(5)=234.95544499d0+1934.136261972222d0*t-0.002075611111d0*t**2      !N'
     >	-0.000000213944d0*t**3+0.000000164972d0*t**4
	astr(6)=282.9373409806d0+1.719457666668d0*t**2+0.000456888889d0*t**2     !ps
     >	-0.000001943279d0*t**3-0.0000000033444d0*t**4
	astr(1)=360.d0*tt-astr(2)+astr(3)
      astr=dmod(astr*RAD,2.d0*pi)
	tdh=0.d0
	do i=1,nn
	  n=nint(fes(i,2));m=nint(fes(i,3))
	  k=n*(n+1)/2+m
	  if(n<1)goto 9005
 	  th1=0.d0;td=nint(fes(i,1));bs=100000
	  do j=1,6
	    dod(j)=td/bs;td=td-bs*dod(j);bs=bs/10
	    if(j>1)dod(j)=dod(j)-5
	  enddo
	  do j=1,6
	    th1=th1+dble(dod(j))*astr(j)
	  enddo
        call BiasTide(fes(i,1),bias)
        call CalcTidefu(mjd,fes(i,1),df,du)
        th1=dmod(th1+(bias+du)*RAD,2.d0*pi)  
        !IERS2010(6.20)
        cp=fes(i,4)*dsin(fes(i,5)*RAD);cn=fes(i,6)*dsin(fes(i,7)*RAD)
        sp=fes(i,4)*dcos(fes(i,5)*RAD);sn=fes(i,6)*dcos(fes(i,7)*RAD)
        !IERS2010(6.19)
        ksp=cp*dcos(th1+dble(m)*rlon)+sp*dsin(th1+dble(m)*rlon)
        ksn=cn*dcos(th1-dble(m)*rlon)+sn*dsin(th1-dble(m)*rlon)
        tdh=tdh+(ksp+ksn)*pnm(k+1)*df
9005	  continue
	enddo
      return
      end
