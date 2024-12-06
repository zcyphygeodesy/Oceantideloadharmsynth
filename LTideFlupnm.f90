       subroutine LTideFlupnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr)
      implicit none
	integer::maxn,nn,n,m,kk,i,j
	real*8::cnm((maxn+2)**2),snm((maxn+2)**2),pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2)
	real*8::gm,tdn(14),tn(14),gr,GRS(6),cosml,sinml,kp
	real*8::BLH(3),rln(3),NFD(5),t,u
 	real*8::rr,rlat,rlon,ae,pi,RAD,flv(4000,3),dk,dh,dl
!---------------------------------------------------------------------------
      kp=1.d0
    	tdn(1:14)=0.d0;gm=GRS(1);ae=GRS(2);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
 	rr=rln(1);rlat=rln(2);rlon=rln(3)
      t=dsin(rlat*RAD);u=dcos(rlat*RAD)
      do n=1,maxn  
         dh=flv(n,1);dl=flv(n,2);dk=flv(n,3);tn=0.d0
         do m=0,n
            j=n*(n+1)/2+m
            cosml=dcos(rlon*RAD*dble(m));sinml=dsin(rlon*RAD*dble(m))
	      tn(1)=tn(1)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*(1.d0+dk)
	      tn(2)=tn(2)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*(kp+(2.d0*dh-dble(n+1)*dk)/dble(n))
	      tn(3)=tn(3)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*(kp+dk)
	      tn(4)=tn(4)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+1)*(1.d0+dk-dh)
	      tn(5)=tn(5)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+1)*(1.d0+dk-dh)
	      tn(6)=tn(6)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+1)*(1.d0+dk)
	      tn(7)=tn(7)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+1)*(1.d0+dk)
	      tn(8)=tn(8)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+1)*dl
	      tn(9)=tn(9)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+1)*dl
	      tn(10)=tn(10)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*dh
	      tn(12)=tn(12)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*(kp+dk)
	      tn(13)=tn(13)+(cnm(j)*cosml+snm(j)*sinml)*dpt2(j+1)*(1.d0+dk)
	      tn(14)=tn(14)+dble(m**2)*(cnm(j)*cosml+snm(j)*sinml)*pnm(j+1)*(1.d0+dk)
         enddo
         tn=tn*dexp(dble(n)*dlog(ae/rr))
 	   tdn(1)=tdn(1)+tn(1)
	   tdn(2)=tdn(2)+tn(2)*dble(n+1)
	   tdn(3)=tdn(3)+tn(3)*dble(n+1)
	   tdn(4)=tdn(4)+tn(4)
	   tdn(5)=tdn(5)+tn(5)
	   tdn(6)=tdn(6)+tn(6)
	   tdn(7)=tdn(7)+tn(7)
	   tdn(8)=tdn(8)+tn(8)
	   tdn(9)=tdn(9)+tn(9)
 	   tdn(10)=tdn(10)+tn(10)
 	   tdn(12)=tdn(12)+tn(12)*dble(n+1)*dble(n+2)
	   tdn(13)=tdn(13)+tn(13)
 	   tdn(14)=tdn(14)+tn(14)
	enddo
	tdn(1)=tdn(1)*gm/rr/gr*1.d3
	tdn(2)=tdn(2)*gm/rr**2*1.0e8
	tdn(3)=tdn(3)*gm/rr**2*1.0e8
	tdn(4)=tdn(4)*gm/rr**2/gr/RAD*36.d5
	tdn(5)=tdn(5)*gm/rr**2/gr/u/RAD*36.d5
	tdn(6)=tdn(6)*gm/rr**2/gr/RAD*36.d5
	tdn(7)=tdn(7)*gm/rr**2/gr/u/RAD*36.d5
	tdn(8)=-tdn(8)*gm/rr/gr/u*1.d3
 	tdn(9)=-tdn(9)*gm/rr/gr*1.d3
	tdn(10)=tdn(10)*gm/rr/gr*1.d3
	tdn(12)=gm/rr**3*tdn(12)*1.d12 !mE
	tdn(13)=-gm/rr**3*tdn(13)*1.d12 !mE
	tdn(14)=gm/rr**3*tdn(14)/u**2*1.d12 !mE
      tdn(11)=tdn(10)-tdn(1)  !地面站点正（常）高变化（mm）
 	return
	end
