      program Otdloalharmsynth
      implicit none
	character*800::outfl,fini
	character*800::line,str
	character*4 fh
	integer i,j,k,n,nn,IY,IM,ID,ih,imn,kk,sn,kln,maxn,astat(5)
      real*8 mjd,sec,tdn(14),cnm(900000),snm(900000),rec(800)
	real*8 GRS(6),flv(4000,3),BLH(3),tmp,pi,RAD,st(7)
	real*8 tm1,tm2,mjd1,mjd2,mjd01,mjd02,tdh
	real*8 plon,plat,phgt,bgntm,endtm,tmdlt,rln(3),NFD(5),gr,t
	integer::status=0
	real*8,allocatable::fes(:,:)
	real*8,allocatable::pnm(:),dpt1(:),dpt2(:)
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
	pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      !Input longitude (degree decimal), latitude (degree decimal), height (m) relative to the sea surface
      !输入经纬度正常高
	plon=121.24d0; plat=29.4281; phgt=17.830
      !Input starting time (long integer time), ending time, time interval (minute) 输入起止时间与时间间隔
      bgntm=20180701;endtm=20180704; tmdlt=30.d0/1440.d0
	maxn=120!Input maximum trucated degree 输入负荷球谐系数最大计算阶数
      !Read FES2004 model 读FES2004格式海潮模型
 	allocate(fes((maxn+2)**2*40,7), stat=astat(1))
	if (sum(astat(1:1)) /= 0) goto 902!内存不足
      open(unit=8,file="FES2004S1.dat",status="old",iostat=status)!Replaceable model here 模型替换
      if(status/=0) goto 902 
      read(8,'(a)') line
      read(8,'(a)') line   !FES2004(11)
      read(8,'(a)') line   !FES2004(11)
      i=0;fes=0.d0
      do while(.not.eof(8))  
	  read(8,*,end=903)st(1),fh,(st(j),j=2,3),tmp,tmp,tmp,tmp,(st(j),j=4,7)
        if(st(2)<maxn+1)then
          i=i+1;st(1)=st(1)*1.d3;fes(i,1:7)=st(1:7);
        endif
	enddo
903   close(8)
      nn=i
	!Read load love numbers 读负荷勒夫数
      flv=0.d0
      open(unit=8,file="Love_load_cm.dat",status="old",iostat=status)
      if(status/=0) goto 902 
      do i=1,6
        read(8,'(a)') line
      enddo
      n=0
	do while(.not.eof(8).and.n<3600)
        n=n+1
	  read(8,*,end=904)i,(flv(n,j),j=1,3)
      enddo
904   close(8)
 	allocate(pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2))
      call tmcnt(bgntm,IY,IM,ID,ih,imn,sec)
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd01=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      call tmcnt(endtm,IY,IM,ID,ih,imn,sec)
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd02=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      if(mjd02<mjd01) goto 902
      open(unit=10,file='reslt.txt',status="replace") !Output file 输出文件
      BLH(1)=plat;BLH(2)=plon;BLH(3)=phgt
      write(10,'(a8,2F12.6,F10.3,F15.6)')'Forcast',plon,plat,phgt,mjd01
      mjd=mjd01;mjd02=mjd02+1.d-6;kk=0
!The Legendre functions are calculated in advance, which can be speeded up.一次性计算勒让德函数，提速
      call BLH_RLAT(GRS,BLH,rln);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      t=dsin(rln(2)*RAD); call BelPnmdt(pnm,dpt1,dpt2,maxn,t)
      do while(mjd<mjd02)
	   call OLoadDFlu(mjd,cnm,snm,maxn,fes,nn)   !The direct influence (cnm,snm) of geopotential coefficient 计算位系数直接影响
         call LTideFlupnm(rln,maxn,cnm,snm,flv,tdn,GRS,pnm,dpt1,dpt2,gr)!Compute ocean tidal load effects 计算海潮负荷效应
         call otidehsynth(mjd,BLH,tdh,fes,nn,maxn,GRS)!Predict sea surface tidal heihgt 预报海洋潮高
         write(10,'(F12.6,20F12.4)')mjd-mjd01,tdh,(tdn(i),i=1,14)!tdh: tidal heihgt (cm), tdh(14): load effects 瞬时潮高(cm)，负荷效应（14列）
         mjd=mjd+tmdlt;kk=kk+1
         if(kk/10*10==kk)write(*, '(a,i9)'), '    Calculated number: ',kk
906      continue
	enddo
      close(10)
      deallocate(fes,dpt1,dpt2)
902	continue
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
!
!************************************************************************
!
      subroutine tmcnt(tm,iyr,imo,idy,ihr,imn,sec)
      !Transform the long integer time (date) agreed by ETideLoad to year, month, day, hour, minute and second
      !ETideLoad格式日期tm转年月日时分秒。
      implicit none
	integer::iyr,imo,idy,ihr,imn,kln
	real*8::tm,sec,tmp,dj0,fd
	character*40::tmstr,astr,aymd,ahms
!-------------------------------------------------------------------------
      write(astr,*)tm
      astr=trim(adjustl(astr));kln=len(astr)
      read(astr(1:4),*)iyr
      read(astr(5:6),*)imo
      read(astr(7:8),*)idy
      ihr=0;imn=0;sec=0.d0
      if(kln>9)read(astr(9:10),*)ihr
      if(kln>11)read(astr(11:12),*)imn
      if(kln>13)read(astr(13:kln),*)sec
      continue
      end

