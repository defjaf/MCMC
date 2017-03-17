ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc TO COMPILE AND LINK WITH ATLAS ccccccccccccccccccccccccccccccccccccccccccccc
cccc g77 -o likico likico.for -L/home/jaffe/local/lib/Linux_Xeon_gcc32_2/ -llapack -lf77blas -lcblas -latlas
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program likelihood_function
      
      call main
      
      end
      
      ! function lnblnk(string)
      ! integer lnblnk, len_trim
      ! character string
      ! lnblnk = len_trim(string)
      ! return
      ! end   
      
      subroutine main
       integer dim,hubblenum,kvalnum,lmax,kmax
       parameter(dim=117,hubblenum=10,kvalnum=7,lmax=10,kmax=40,mmax=41)
       double precision hubble(hubblenum)
       integer kwav(kvalnum),mult(kvalnum)    
       double precision dlnk(hubblenum,kvalnum,lmax-1),
     *apowers(hubblenum,kvalnum,lmax-1),
     *transferf(hubblenum,kvalnum,lmax-1)
       double complex xi(kvalnum,mmax,(kmax+1)**2),alm(dim)
       double precision alikelihood,hub
       common /input/kwav,mult,hubble,dlnk,apowers,transferf,xi,alm

       double precision hmin, hmax
       integer nstep

      character*150 datdir, almfile
      
      datdir = './dat/'
      almfile = ''

      call readdata(datdir, almfile)
c      do i=1,117
c         print*,i,alm(i)
c      enddo
c       do i=1,10
c          print*,hubble(i)
c       enddo
c      print*,dlnk(1,1,1),apowers(1,1,1),transferf(1,1,1)
c      print*,dlnk(10,7,9),apowers(10,7,9),transferf(10,7,9)
c      print*,xi(1,1,1),xi(1,13,169)
c      print*,xi(2,1,1),xi(2,21,441)
c      print*,xi(3,1,1),xi(3,25,625)
c      print*,xi(4,1,1),xi(4,31,961)
c      print*,xi(5,1,1),xi(5,33,1089)
c      print*,xi(6,1,1),xi(6,37,1369)
c      print*,xi(7,1,1),xi(7,41,1681)
c      alpha=0.364298
c      beta=0.1290
c      gamma=1.84952
c      do i=1,21
c      hub=(70.d0-52.d0)/dble(20)*(i-1)+52.d0
c      print*,hub,alikelihood(1.d6,hub,alpha,beta,gamma)
c      enddo 
      alpha=0.
      beta=0.
      gamma=0.
      hmax = 70
      hmin = 52
      nstep = 1+(hmax - hmin)/2


      do i=1,nstep
       hub=(hmax-hmin)/dble(nstep-1)*(i-1)+hmin
c      alpha=2.*acos(-1.)/20*(i-1)
      print*,alpha,hub,alikelihood(1.d0,hub,alpha,beta,gamma)
      enddo

      return
      end

      subroutine readdata(datdir, almfile)
      character*150 datdir, almfile
Cf2py intent(in) datdir, almfile

      integer dim,hubblenum,kvalnum,lmax,kmax
      parameter(dim=117,hubblenum=10,kvalnum=7,lmax=10,kmax=40,mmax=41)
      double precision hubble(hubblenum)
      integer kwav(kvalnum),mult(kvalnum)
      double precision dreal1,dreal2,dlnk(hubblenum,kvalnum,lmax-1),
     *apowers(hubblenum,kvalnum,lmax-1),
     *transferf(hubblenum,kvalnum,lmax-1)
      double complex xi(kvalnum,mmax,(kmax+1)**2),alm(dim)
      common /input/kwav,mult,hubble,dlnk,apowers,transferf,xi,alm

      !!integer len_trim
      !!external len_trim

      if (len_trim(almfile).eq.0) then
          almfile = 'alm64_1.dat'
      endif

      open(22,file=datdir(1:len_trim(datdir))//almfile)
      do i=1,dim
         read(22,*)idum,dreal1,dreal2
         alm(i)=complex(dreal1,dreal2)
      enddo
      close(22)

      open(23,file=datdir(1:len_trim(datdir))//'hubblelist.dat')
      do i=1,hubblenum
         read(23,*)int1
         hubble(i)=dble(int1)
      enddo
      close(23)
      open(25,file=datdir(1:len_trim(datdir))//'transf.dat')
      do i=1,hubblenum
         do j=1,kvalnum
            do ll=1,lmax-1
               read(25,*)idum,idum,dlnk(i,j,ll),apowers(i,j,ll),
     *transferf(i,j,ll)
            enddo
         enddo
      enddo
      close(25)
      open(27,file=datdir(1:len_trim(datdir))//'kvalues.dat')
      do i=1,kvalnum
         read(27,*)kwav(i),mult(i)
c      print*,kwav(i),mult(i)
      enddo
      close(27)
      do i=1,kvalnum
         do j=1,mult(kvalnum)
            do k=1,kwav(kvalnum)+1
               xi(i,j,k)=complex(0.,0.)
            enddo
         enddo
      enddo
      open(28,file=datdir(1:len_trim(datdir))//'xiall.dat')
      do i=1,kvalnum
         do j=1,mult(i)
            do k=1,(kwav(i)+1)**2
            read(28,*)idum,idum,idum,dreal1,dreal2
            xi(i,j,k)=complex(dreal1,dreal2)
            enddo
         enddo
       enddo
       close(28)
       return
       end

       
       function alikelihood(ampl,h0,alpha,beta,gamma)
Cf2py intent(in) ampl, h0, alpha, beta, gamma
       
       complex wigner2
       double complex sumxi,sumlik,almr(117),clm(117,117)
       double precision dlogtemp,alogsum,alikelihood,dlogdet,alik1
       integer ipiv(117)
       double complex work(9360)
       double precision transfer1(10),transfer2(10),tran1,tran2,dtran1,
     *dtran2,h0,dlogampl,ampl
       integer dim,hubblenum,kvalnum,lmax,kmax
       parameter(dim=117,hubblenum=10,kvalnum=7,lmax=10,kmax=40,mmax=41)
       double precision hubble(hubblenum)
       integer kwav(kvalnum),mult(kvalnum)    
       double precision dlnk(hubblenum,kvalnum,lmax-1),
     *apowers(hubblenum,kvalnum,lmax-1),
     *transferf(hubblenum,kvalnum,lmax-1)
       double complex xi(kvalnum,mmax,(kmax+1)**2),alm(dim)
       logical fail
       external wigner2
       common /input/kwav,mult,hubble,dlnk,apowers,transferf,xi,alm

       fail = .false.
       dlogampl=dlog(ampl)
       dlogtemp=dlog(2.726**2*2.d12)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc compute correlation matrix  ccccccccccccccccccccccccccccccccccccccccccc
       do i=1,(lmax+1)**2-4
          do j=1,(lmax+1)**2-4
             clm(i,j)=complex(0.d0,0.d0)
          enddo
       enddo
       do kk=1,kvalnum
       do ll1=2,lmax
          do mm1=-ll1,ll1,1
             ind1=ll1**2+ll1+mm1-3
             do ll2=2,lmax
                if(kwav(kk).ge.ll1.and.kwav(kk).ge.ll2)then
                   do i=1,hubblenum
c                  print*,hubble(i)
                   transfer1(i)=transferf(i,kk,ll1-1)
                   transfer2(i)=transferf(i,kk,ll2-1)
                   enddo
                 call polint(hubble,transfer1,hubblenum,h0,tran1,
     *                       dtran1,fail)
                 if (fail) then
                    print *, 'Failure in polint (1); h0=',h0
                    alikelihood = 0.0
                    return
                 endif
                 call polint(hubble,transfer2,hubblenum,h0,tran2,
     *                       dtran2,fail)
                 if (fail) then
                    print *, 'Failure in polint (2); h0=',h0
                    alikelihood = 0.0
                    return
                 endif
c                 print*,kwav(kk),ll1,ll2,tran1,tran2
                endif
                do mm2=-ll2,ll2,1
                   ind2=ll2**2+ll2+mm2-3
                   sumxi=complex(0.d0,0.d0)
                   do i=1,mult(kk)
                      sumxi=sumxi+xi(kk,i,ind1+4)*conjg(xi(kk,i,ind2+4))
                   enddo
c                   print*,ind1,ind2,sumxi
                   if(kwav(kk).ge.ll1.and.kwav(kk).ge.ll2)then
                   alogsum=dlog(dlnk(1,kk,ll1-1))+
     *                     dlog(apowers(1,kk,ll1-1))+dlogtemp+dlogampl
                   clm(ind1,ind2)=clm(ind1,ind2)+exp(alogsum)*
     *                            tran1*tran2*sumxi*dble(120)
                   endif
                enddo
             enddo
          enddo
       enddo
       enddo
c       open(12,file=datdir(1:len_trim(datdir))//'clmc.dat')
c       do i=1,117
c       do j=1,117
c       write(12,*)i,j,dreal(clm(i,j)),dimag(clm(i,j))
c       enddo
c       enddo
c       close(12)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc invert correlation matrix  cccccccccccccccccccccccccccccccccccccccccc
       call zgetrf(117,117,clm,117,ipiv,info1)
c       print*,info1
       dlogdet=0.d0                                     !zgetrf returns as clm the LU decomposition
       do i=1,117                                       !of clm. Therefore the absolute value of the
          dlogdet=dlogdet+dlog(abs(clm(i,i)))           !determinant is the product of the absolute values of
       enddo                                            !the diagonal of the matrix containing the LU factors
       dlogdet=dlogdet+dlog(2.d0*dacos(-1.d0))
       call zgetri(117,clm,117,ipiv,work,9360,info)
c       print*,info
c       open(13,file=datdir(1:len_trim(datdir))//'invclmc.dat')
c       do i=1,117
c       do j=1,117
c       write(13,*)i,j,dreal(clm(i,j)),dimag(clm(i,j))
c       enddo
c       enddo
c       close(13)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc rotate alms  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       do i=1,dim
          almr(i)=complex(0.d0,0.d0)
c          print*,i,alm(i)
       enddo
       do ll=2,10
          ind1=ll**2-3
          ind2=ll**2+2*ll-3
          do mm=-ll,ll,1
             ind=ll**2+ll+mm-3
             ii=0
             do i=ind1,ind2
                almr(ind)=almr(ind)+
     *          wigner2(alpha,beta,gamma,ll,-ll+ii,mm)*alm(i)
                ii=ii+1
c                print*,ind,almr(ind)
             enddo
c          print*,ind,almr(ind)
          enddo
       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc  compute -1/2*conjugate(alm)*inv(C_lml'm')*alm ccccccccccccccccccccccccc
       sumlik=complex(0.d0,0.d0)
       do i=1,117
          do j=1,117
             sumlik=sumlik+conjg(almr(j))*clm(j,i)*almr(i)
          enddo
       enddo
c       print*,sumlik
       alik1=-1.d0/2.d0*dreal(sumlik)          !the exponent
       alik1=alik1-1.d0/2.d0*dlogdet           !the exponent plus the logarithm of 1/sqrt(2pi*det[C_lml'm'])
       alikelihood=alik1                       !the logarithm of the likelihood
c       print*,flik
       return
       end


      complex function wigner2(alpha,beta,gamma,j,m1,m2)
      real alpha,beta,gamma
      integer j,m1,m2
      complex wigner
      external wigner
      wigner2=float((-1)**(m2-m1))*
     *conjg(wigner(alpha,beta,gamma,j,m1,m2))
      return
      end

      complex function wigner(alpha,beta,gamma,j,m1,m2)
      real alpha,beta,gamma
      integer j,m1,m2
      complex epsilon1,epsilon2
      real jacobip
      external jacobip

      mu=abs(m1-m2)
      nu=abs(m1+m2)
      s=float(j)-1./2.*(mu+nu)
c      print*,'mu=',mu,'nu=',nu,'s=',s
      if(m2.ge.m1)then
      ksi=1
      else
      ksi=(-1)**(m2-m1)
c      print*,'ksi=',ksi
      endif
      factorial=1.
      do i=1,nu
         factorial=factorial*sqrt((s+mu+i)/(s+i))
      enddo
c      print*,'f=',factorial
      alittled=ksi*factorial*(sin(beta/2.))**mu*(cos(beta/2.))**nu*
     *        jacobip(int(s),float(mu),float(nu),cos(beta))
      epsilon1=complex(cos(m1*alpha),-sin(m1*alpha))
      epsilon2=complex(cos(m2*gamma),-sin(m2*gamma))
      wigner=epsilon1*alittled*epsilon2
      return
      end

      function jacobip(n,a,b,x)
      real jacobip,a,b,x,jp0,jp1,jp2
      jp0=1.
      jp1=1./2.*(2.*(a+1.)+(a+b+2)*(x-1.))
      jp2=1./8.*(4.*(a+1.)*(a+2.)+4.*(a+b+3.)*(a+2.)*(x-1.)+
     *    (a+b+3.)*(a+b+4.)*(x-1.)**2)
      if(n.le.2)then
       if(n.eq.0)then
       jacobip=jp0
       endif
       if(n.eq.1)then
       jacobip=jp1
       endif
       if(n.eq.2)then
       jacobip=jp2
       endif
      else
      pj1=jp1
      pj2=jp2
      do l=3,n
         pj=((a+b+float(2*l)-1.)*(a**2-b**2+x*(a+b+float(2*l)-2.)
     *       *(a+b+float(2*l)))*pj2
     *      -2.*(float(l)+a-1.)*(float(l)+b-1.)*(float(2*l)+a+b)*pj1)/
     *       (2.*float(l)*(a+b+float(l))*(a+b+float(2*l)-2))
         pj1=pj2
         pj2=pj
      enddo
      jacobip=pj
      endif
      return
      end


      subroutine polint(xa,ya,n,x,y,dy,fail)
      integer n,NMAX,l
      double precision dy,x,y,xa(n),ya(n)
      parameter(NMAX=20)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      logical fail
      fail=.false.
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
            dift=abs(x-xa(i))
            if(dift.lt.dif)then
               ns=1
               dif=dift
            endif
            c(i)=ya(i)
            d(i)=ya(i)
 11   enddo
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
            do 12 i=1,n-m
                  ho=xa(i)-x
                  hp=xa(i+m)-x
                  w=c(i+1)-d(i)
                  den=ho-hp
                  if(den.eq.0.d0)then
                     print*, 'failure in polint: n, x, y, dy=',n,x,y,dy
                     fail=.true.
                     return
                  endif
                  den=w/den
                  d(i)=hp*den
                  c(i)=ho*den
 12         enddo
            if(2*ns.lt.n-m)then
               dy=c(ns+1)
            else
               dy=d(ns)
               ns=ns-1
            endif
            y=y+dy
 13   enddo
      return
      end


      SUBROUTINE ratint(xa,ya,n,x,y,dy,fail) 
      INTEGER n,NMAX 
      double precision dy,x,y,xa(n),ya(n),TINY 
      PARAMETER (NMAX=10,TINY=1.d-25)
      INTEGER i,m,ns 
      REAL dd,h,hh,t,w,c(NMAX),d(NMAX) 
      logical fail
      fail=.false.
      ns=1 
      hh=abs(x-xa(1)) 
      do 11 i=1,n 
            h=abs(x-xa(i)) 
            if (h.eq.0.d0)then 
                y=ya(i) 
                dy=0.d0 
                return 
            else if (h.lt.hh) then 
                     ns=i 
                     hh=h 
            endif 
            c(i)=ya(i) 
            d(i)=ya(i)+TINY
 11    enddo
       y=ya(ns) 
       ns=ns-1 
       do 13 m=1,n-1 
          do 12 i=1,n-m
                w=c(i+1)-d(i) 
                h=xa(i+m)-x
                t=(xa(i)-x)*d(i)/h
                dd=t-c(i+1)
                if(dd.eq.0.d0) then
                    print *,'failure in ratint: n, x, y, dy=',n,x,y,dy 
                    fail = .true.
                    return
                endif
                dd=w/dd 
                d(i)=c(i+1)*dd 
                c(i)=t*dd 
 12        enddo 
           if (2*ns.lt.n-m)then 
               dy=c(ns+1) 
           else 
               dy=d(ns) 
               ns=ns-1 
           endif 
       y=y+dy 
 13    enddo  
       return 
       END
