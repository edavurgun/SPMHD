PROGRAM AxisSPHYNX_MHD

      implicit none

      integer, parameter :: n1=1500000
      integer, parameter :: n2=2*n1
      integer, parameter :: nnl=1000001
      integer, parameter :: nBhmean=100
      integer, parameter :: ni0r=60
      integer, parameter :: vecMax=150
      integer, parameter :: nt=3*n2!2*n2
      integer, parameter :: ntree=2*n2!n2
      integer, parameter :: iodelay=500
      integer, parameter :: initial_adjust_iterations=10
      integer, parameter :: nwrite=0
      integer, parameter :: itmaxNR=6   !maximum iterations in Newton-Raphson (temperature)
      logical, parameter :: averVelr=.false.
      logical, parameter :: crossedMom=.true.
      logical, parameter :: timeaxis=.false.
      logical, parameter :: timeTemp=.true.
      logical, parameter :: balsara_lim=.true.
      logical, parameter :: cleandiv=.true.
      logical, parameter :: derivphi1=.true.    !picks the SPH flavor of derivatives. 
      !For PSR-Disk problem .true. gives better deriv_vphi_z (but -r-derivatives are similar and not very good)
      logical, parameter :: profiling=.false.

      double precision, parameter :: pi=acos(-1.d0)
      double precision, parameter :: pim1=1.d0/pi
      double precision, parameter :: pi2=pi*0.5d0
      double precision, parameter :: year=365.d0*86400.d0
      double precision, parameter :: pc=3.086d18      ! parsec to cm
      double precision, parameter :: g=6.6743d-08
      double precision, parameter :: clight=3.0d10
      double precision, parameter :: P_PSR=1.69d-03
      double precision, parameter :: Bp=9.6d7
      double precision, parameter :: phi0=pi/2.d0
      double precision, parameter :: psi=(0.d0*pi/180.d0)
      double precision, parameter :: msun=1.989d33
      double precision, parameter :: rhoCrit=2.d-12
      double precision, parameter :: TCrit=2.d3
      double precision, parameter :: rns=1.5d6
      double precision, parameter :: masaPSR=1.4*msun
      double precision, parameter :: rgas=8.317d7
      double precision, parameter :: arad=7.564d-15
      double precision, parameter :: mui=0.615
      double precision, parameter :: rs=2.5d6
      double precision, parameter :: xlim=2.d6     !critical distance from axis 
      double precision, parameter :: ylimpos=5.d6     !critical distance up 
      double precision, parameter :: ylimneg=-5.d6     !critical distance  down
!        to remove particles from the simulation box
      double precision, parameter :: gamma=5.d0/3.d0     ! for the barotropic EOS, internal energy
      double precision, parameter :: gammam1=gamma-1.d0
      !RESTART parameters
      integer, parameter :: lstart=1
      double precision, parameter :: initial_timestep=1.d-8
      double precision, parameter :: timeinit=0.d0
      character*21, parameter     :: inputfile='D05LC.07999'
      !End RESTART parameters

      double precision, parameter :: Cour=0.4d0
      double precision, parameter :: C_accel=0.25d0
      double precision, parameter :: tol_temp=0.05d0
      double precision, parameter :: tol_taxis=0.05
      double precision, parameter :: alfaAV=1.d0
      double precision, parameter :: betaAV=2.d0
      double precision, parameter :: alfau=0.05
!      double precision, parameter :: alfau=0.5d0
      double precision, parameter :: alfaB=1.d0
      double precision, parameter :: hindex=5.d0
      double precision, parameter :: epsdiv=0.0d0  !idefault 0.2
      double precision, parameter :: mass_scaling=1.d0
      double precision, parameter :: mix=0.5d0   !extent of mixing in tensile instability (0,1) mix=-1 also posible
      double precision, parameter :: rboundary_scaling=1.d0     ! fine tuning of bounday condition
      double precision, parameter :: minvalue=1.d-60
      double precision, parameter :: tolNR=1.d-05
      double precision, parameter :: rpos=1.d7
      double precision, parameter :: zpos=0.d0
      double precision, parameter :: radpart=0.7d5
      double precision, parameter :: rpos2=8.0d7
      double precision, parameter :: radpart2=1.d5
      double precision, parameter :: rpos3=1.5d8
      double precision, parameter :: radpart3=1.d5

      double precision, parameter :: mu0=4.d0*pi    ! Magnetic permeability. In IS: 4*pi*1e-7. In cgs:  4*pi (???)
!      double precision, parameter :: mu0=1.d0       ! Magnetic permeability in code units. B(cgs)= B(code) x sqrt(4 pi)
      double precision, parameter :: fclean=1.d0    ! To clean divergence. Better <=1
      double precision, parameter :: sigma_c=1.d0   ! Adimensional constant to control the magnetic  decay term during cleaning

      logical remove
      integer signature,value1,signatureC,signatureNC,index,irhomax,irhomin
      integer nvmin,nvmax,nvr,nsatura,nje,nm,npp,nnm,i1,i2,npp0,nmt
      integer invmax,nvSum,npointer1,ntotal,nmmm,ncount,nodo,ncount1,ncount2,ncount3
      integer nnma,l,kl,i,invmin,j,jg,jVec,k,n,i_out,n11,n22

      integer, dimension(n1) :: nv,nvvv
      integer, dimension(n2) :: ngn,ngh
      integer, dimension(nt) :: nnod,npa,np,npn,npart,npointer
      integer, dimension(0:nt) :: nh1,nh2
      integer, dimension(ntree) :: naux,nauxa
      integer, dimension(n1,vecMax) :: vecMatrix
      integer, dimension(6)    :: part,part2,part3

      double precision m1,nu2,nindex1,nindex2,pkv,pmt,ro3dij,ro3dmax,ro3dmin
      double precision momentumxx,momentumyy,momentumzz,momentumxxp,momentumyyp,romax
      double precision momentumxxd,momentumyyd,momentumzzd,momentumxxl,momentumyyl,momentumzzl
      double precision momentumxxa,momentumyya,momentumzza

      double precision kernxxi,kernxyi,kernyxi,kernyyi,ux,uy,v1,val,rhomax,rr2,rr2_32,rr,rrr
      double precision kernxxj,kernxyj,kernyxj,kernyyj,kernxx,kernyy,kx
      double precision accel,angmomphi,angmomtot,angmomr,angmomz
      double precision lambda,masstot,lcentered,ax,ay,Bcx,Bcz,Bcphi,Bmo12,Bmod
      double precision Brad,Btheta,Bphi,Bradav,Bthetaav,Bphiav
      double precision dtnma,dtnme,tt,dt,tt1,xt,dx,dy,h1,d02,d05,pi2v1,w1d,dter
      double precision Bmod2,txx,txy,tyy,hhi,d11,d21,rij_2,rdenom,rinvdenom
      double precision ddmax,ddist,rad,h0,h11,hexp,hfac
      double precision divi,divj,roti,rotj,fbalsi,fbalsj,fij,qiji,qijj,qijiu,qijju,tqij
      double precision proi,proj,proeneri,tp,ti3,fix0,djix,djiy,vjix,vjiy,vjiz
      double precision gradxha,gradyha,gradxhb,gradyhb,checkdeltar
      double precision dteri,dvx,dvy,dmy,dmy1,dmy2,smdiv,tau_a,beta,hhj,h2,uijx,uijy
      double precision v2,pi2v2,w2d,dterj,dj,vijrij,r3dij,rij,wij,vijsignal,vijsignalu
      double precision dmy01,dmy02,dmy11,dmy12,cm,cm0,d1,d2,dc,di,dtn,du,gradxhab,gradyhab
      double precision aviscox,aviscoy,dpx,dpy,dpxavu,dpyavu,Bij2,Fabha,Fabhb
      double precision vijsignalBi,vijsignalBj,dmyr,dmyz,dmyphy,vijsignaldis
      double precision viscBdis,viscBi,viscBj,dmyB,dmyBclean,Bdist,Btot2,hoop_r,hoop_phi
      double precision vijsignalab,resistiv_ab,Jdis
      double precision hoop,ratio,a11,a12,a13,a21,a22,a23,a31,a32,a33
      double precision smdi,egrav,eint,ekin,elost,emag,errdivB,errdivBav,ycma
      double precision time_accel,time_temp,time_Cour,time_axis,errdivBmax,errdivBmin,etot0
      double precision totalmassC,space,rdmax,zdmaxneg,zdmaxpos,rtr
      double precision rtzu,rtzd,aux1,aux2,aux3,dpdro,desplr,desplz,rtr5,rtz5,xcm,ycm,dmass
      double precision Sij_i11,Sij_i12,Sij_i13,Sij_i21,Sij_i22,Sij_i23,Sij_i31,Sij_i32,Sij_i33
      double precision Sij_j11,Sij_j12,Sij_j13,Sij_j21,Sij_j22,Sij_j23,Sij_j31,Sij_j32,Sij_j33
      double precision vij1,vij2,vij3,vmod,vr,x1,xcma,y1,ycm_check,angmomtot_0
      double precision dudt1,dtemp,gnr,dgnr,ueos
      double precision rLC,omegaPSR,dist,dist2,dist3,polar_ang,xdisk,dBr,dBtheta,dBphi,deltat,Br,Bz
      double precision tempaver,densaver,tempaver2,densaver2,tempaver3,densaver3
      double precision tempaver_0,tempaver2_0,tempaver3_0,densaver_0,densaver2_0,densaver3_0
      
      double precision, dimension(4)    :: ccf0
      double precision, dimension(n1)   :: sum3,evisco,eviscu,maxvsignal
      double precision, dimension(n1)   :: divB,cha0,cha
      double precision, dimension(n1)   :: checkInorm,checkdeltax,checkdeltay
      double precision, dimension(n1)   :: divvphir,divvphiz,divvrr,divvrz,divvzr,divvzz
      double precision, dimension(n1)   :: gradrBphi,gradrBr,gradrBz
      double precision, dimension(n1)   :: proener,sum1B,sum2B,sum3B,B2mu0,phi,eviscoB,eviscoBbis
      double precision, dimension(n1)   :: ap1,ap2,ap1_2,ap2_2,ri_2,ap1p2,ap2p2
      double precision, dimension(n1)   :: lz,lz0,wz
      double precision, dimension(n1)   :: avevel,bfieldr,gradPr
      double precision, dimension(n1)   :: fgravr,fgravz
      double precision, dimension(n1)   :: Jdis_r,Jdis_z,Jdis_phi


      double precision, dimension(ionmax) :: am,zm
      double precision, dimension(n2,2) :: a,a0,f
      double precision, dimension(n1,3) :: Bdivclean,Bindis
      double precision, dimension(n2,3) :: v0,v,B, B0


      double precision, dimension(n2)   :: sum1,sum2,sum5,va_sroot,ro2d,ro3d,h,mass,massa,masaa
      double precision, dimension(n2)   :: cxx,cxy,cyy,c,p,u,u0,temp,temp0,pro,sig,x,y
      double precision, dimension(n2)   :: divv,curlz,vol,valfven2,phidiv,phidiv0,dpdt,dudt,dudv
      double precision, dimension(n1)   :: abar,zbar
      double precision, dimension(nt)   :: dxx,xcet,ycet,xce,yce

      double precision, dimension(ntree,1)   :: xia,yia
      double precision, dimension(ntree,2)   :: xi,xd,yi,yd

      double precision, dimension(n2,3,3)  :: Sij
      logical, dimension(nt) :: particle
      logical, dimension(n2) :: axis
      logical, dimension(n1) :: ready
  
      character*14 nomfs1,nomfs2
      character*9  prefix
      character*5  suffix

!        ---------------------------------------------------------
!               redefining total number of particles so that 
!               we can remove particles from the simulation box
!               as the calculation progresses.
!       ---------------------------------------------------------
                 n11=n1
                 n22=n2

!     Interpolating kernels are the 2d-sinc kernels with index n (default n=5)
!     Axis corrections via negative ghosts particles
      data (ccf0(i),i=1,4)/1.3090244919d-01,1.935848488354d-02,-6.164290621127d-03,5.22450269233d-02/

      pkv=ccf0(1)*hindex+ccf0(2)/hindex+ccf0(3)/hindex**2+ccf0(4)
      dtnma=initial_timestep
      dtnme=dtnma
      dtn=(dtnma+dtnme)*0.5d0
      tt=timeinit

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      open(1,file='D05LC.07999')
      open(2,file='tree_datax05B8T00_00')
      open(3,file='time_stepx05B8T00_00')
      open(9,file='Particle1x05B8T00_00.dat')

      v  = 0.d0
      phidiv  = 0.d0
      phidiv0 = phidiv
      B=0.d0
      sum1 = 0.d0
      sum2 = 0.d0
      sum3 = 0.d0

      ncount=0
      ncount1=0
      ncount2=0
      ncount3=0

      if(profiling) call profile(0,'')
      !read data
       do i=1,n11
        read(1,*) a(i,1),a(i,2),rrr,mass(i),h(i),&
          &    c(i),temp(i),ro3d(i),ro2d(i),v(i,1),&
          &    v(i,2),v(i,3),u(i),p(i),wz(i),sum1(i),sum2(i),fgravr(i),&
          &    fgravz(i),dmy,dmy,checkInorm(i),&
           &    nv(i)
           dist=sqrt((a(i,1)-rpos)**2+(a(i,2)-zpos)**2)
           if(dist.le.radpart .and. ncount1.le.5) then
              ncount1=ncount1+1
              part(ncount1)=i
           endif  
            dist2=sqrt((a(i,1)-rpos2)**2+(a(i,2)-zpos)**2)
           if(dist2.le.radpart2 .and. ncount2.le.5) then
              ncount2=ncount2+1
              part2(ncount2)=i
            !  print*, dist2,ncount2,i
           endif
           dist3=sqrt((a(i,1)-rpos3)**2+(a(i,2)-zpos)**2)
           if(dist3.le.radpart3 .and. ncount3.le.5) then
              ncount3=ncount3+1
              part3(ncount3)=i              
       !       print*, dist3,ncount3,i
           endif
       enddo
       if(profiling) call profile(1,'   read_data')
       ! Assign initial values
         u0=u
         temp0=temp
         v0=v 
        do i=1,n11                 ! specific angular momentum of the disc
         v(i,3)=wz(i)*a(i,1)
         v0(i,3)=v(i,3)
         lz0(i)=v(i,3)*a(i,1)
         lz(i)=lz0(i)
        enddo  

       ! Set initial values of the magnetic field. Determinationm of 
       ! the harmonic mean of the three components of the magnetic field
       ! which is the one the disk 'sees' (? Eksi~Alpar)after several pulsar periods
         omegaPSR=2.d0*pi/P_PSR
         rLC=clight/omegaPSR
 !!         go to 4576
!!!!       !$omp parallel private(i,j,Bradav,Bthetaav,Bphiav,dist,xdisk,polar_ang,&
!!!!       !$omp phi,Brad,Btheta,Bphi)
!!!!       !$omp do 
        do i=1,n11 
          Bradav=0.d0
           Bthetaav=0.d0
            Bphiav=0.d0
          dist=dsqrt(a(i,1)**2+a(i,2)**2)
          xdisk=dist/rLC 
          if(a(i,2).ne.0.d0) then
          polar_ang=atan(abs(a(i,1)/a(i,2)))
         else
             polar_ang=pi/2.d0
         endif
         if(a(i,2).lt.0.d0) polar_ang=pi-polar_ang
           deltat=P_PSR/nBhmean
         do j=1,nBhmean+1
           tt1=(j-1)*deltat
        call  Bdeutsch(psi,xdisk,polar_ang,rLC,Bp,P_PSR,rns,&
              &   omegaPSR,phi0,tt1,Brad,Btheta,Bphi)
!          But we have to pass from spherical (r,theta) to axisymmetric coordinates (ro,z)
           Br=Brad*dsin(polar_ang)+ Btheta*dcos(polar_ang)
           Bz=Brad*dcos(polar_ang)-Btheta*dsin(polar_ang)
           Bphi=Bphi
           Bradav=Bradav+Br*Br
           Bthetaav=Bthetaav+Bz*Bz
           Bphiav=Bphiav+Bphi*Bphi
          enddo
!            Note that magnetic field is enforced to have positive components
            B(i,1)=sqrt(Bradav/(nBhmean+1))
            B(i,2)=sqrt(Bthetaav/(nBhmean+1))
            B(i,3)=sqrt(Bphiav/(nBhmean+1))
         enddo
!!!!         !$omp end do 
!!!!         !$omp end parallel
!! 4576     continue
       do i=1,n11
         B0(i,1) = B(i,1)
         B0(i,2) = B(i,2)
         B0(i,3) = B(i,3)
       enddo

       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       !Maximum distance in r (radial)  and z (vertical)
       rdmax=maxval(abs(a(:,1)))
       zdmaxpos=maxval(a(:,2))
       zdmaxneg=abs(minval(a(:,2)))

       rtr=rdmax
       rtzu=abs(zdmaxpos)
       rtzd=abs(zdmaxneg)

    
       !Ghosts particles
       desplr=0.d0
       desplz=0.d0
       rtr5=rtr/5.d0         !Copy just a fifth of the domain to save memory
       rtz5=rtzu/5.d0

       do i=1,n22
         ngh(i)=i
         axis(i)=.false.
       enddo

       do i=1,nc
         call ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,i,n,&
                  &  a,v,B,n11,n22,axis,ngh)
       enddo
       desplr=1.5d9 
       desplz=1.5d9
       a(:,1)=a(:,1)+desplr
       a(:,2)=a(:,2)+desplz
       a0=a
       v0=v
       do i=1,n
         h(i)=h(ngh(i))
         u(i)=u(ngh(i))
         mass(i)=mass(ngh(i))
         ro3d(i)=ro3d(ngh(i))
       enddo
       u0=u
       !*** CENTER OF MASS ***
       ycm=0.d0
       dmass=0.d0
       do i=1,n11
         ycm=ycm+mass(i)*a(i,2)
         dmass=dmass+mass(i)
       enddo
       !ycm = ycm/dmass
       xcm=desplr
       ycm=desplz

       !*** SIZE OF THE SYSTEM ***
       ddmax = 0.0d0
       do i = 1,n
         ddist=sqrt((xcm-a(i,1))**2+(ycm-a(i,2))**2)
         ddmax=max(ddmax,ddist)
       end do
       rad = ddmax

!     *** Pressure (ideal gas) + Radiation ***
         do  i=1,n11
             p(i)=rgas*ro3d(i)*temp(i)/mui+arad/3.d0*temp(i)**4
             u(i)=1.5*rgas*temp(i)/mui+arad*temp(i)**4/ro3d(i)
             c(i)=dsqrt(rgas*temp(i)/mui)
             u0(i)=u(i)
        enddo
          do i=1,n22   !iInitializing Particle-ghosts translator
            ngh(i)=i
          axis(i)=.false.
           enddo

       !     **************************************************************
       !     Iterations start, nnl is the total number of iterations
       !     kl=1,2 is centering the schemei
       do l=lstart,nnl
         do kl=1,2
           dt = (dtnma / 2.d0) * dble(kl)
           if(kl.eq.1) write(3,'(1x,i5,7(1x,e14.7))') l,time_accel,time_temp,time_Cour,time_axis,dtnma,tt
           ! TREE CALCULATION
           if(profiling) call profile(0,'')
           do i = 1,n
             massa(i)=mass(i)
             x(i)=a(i,1)
             y(i)=a(i,2)
             npointer(i)=i
           end do
           do i=1,4
             nauxa(i)=0
           end do
           nnma=0
           ntotal=0
           npn(1)=1
           npn(2)=n
!           xt=pi*sqrt(xcm**2+ycm**2)/sqrt(7.d0)
            rdmax=maxval(abs(a(:,1)))    
            xt=2.*(sqrt(2.d0)-0.3d0)*rdmax     !dynamic max-size of the box
           nje=0
           nm=1
           xia(1,1)=0.d0
           yia(1,1)=0.d0
     6     nje=nje+1
           npp=0
           dx=xt/2.d0**nje
           dy=dx
           do j=1,nm
             np(2*j-1)=npn(2*j-1)
             np(2*j)=npn(2*j)
             xi(j,1)=xia(j,1)
             yi(j,1)=yia(j,1)
             naux(j)=nauxa(j)
           enddo
           nnm=nm
           if(nnm .ge. nnma) nnma=nnm
           if(nnm .gt. ntree-10) nsatura=1
           nm=0
           do j=1,nnm
             nmmm=0
             do i1=1,2
               xi(j,i1)=xi(j,1)+(i1-1)*dx
               xd(j,i1)=xi(j,i1)+dx
               do i2=1,2
                 yi(j,i2)=yi(j,1)+(i2-1)*dy
                 yd(j,i2)=yi(j,i2)+dy
                 npp0=0
                 xcma=0.d0
                 ycma=0.d0
                 pmt=0.d0
                 do i=np(2*j-1),np(2*j)
                   if((xi(j,i1) .le. x(i)) .and. (x(i).lt.xd(j,i1))) then
                     if((yi(j,i2) .le. y(i)) .and. (y(i).lt.yd(j,i2))) then
                       npp=npp+1
                       npp0=npp0+1
                       xcma=xcma+x(i)*massa(i)
                       ycma=ycma+y(i)*massa(i)
                       pmt=pmt+massa(i)
                       x1=x(npp)
                       y1=y(npp)
                       m1=massa(npp)
                       npointer1=npointer(npp)
                       x(npp)=x(i)
                       y(npp)=y(i)
                       massa(npp)=massa(i)
                       npointer(npp)=npointer(i)
                       x(i)=x1
                       y(i)=y1
                       massa(i)=m1
                       npointer(i)=npointer1
                     endif
                   endif
                 enddo
                 if(npp0-1.lt.0) goto 25
                 if(npp0-1.eq.0) goto 26
        27	     ntotal=ntotal+1
                 if(ntotal .gt. nt-50) nsatura=1
                 nm=nm+1
                 nmmm=nmmm+1
                 xcet(ntotal)=(xi(j,i1)+xd(j,i1))/2.d0
                 ycet(ntotal)=(yi(j,i2)+yd(j,i2))/2.d0
                 npa(ntotal)=naux(j)
                 nauxa(nm)=ntotal
                 dxx(ntotal)=dx
                 particle(ntotal)=.false.
                 npn(2*nm-1)=npp-npp0+1
                 npn(2*nm)=npp
                 xia(nm,1)=xi(j,i1)
                 yia(nm,1)=yi(j,i2)
                 goto 28
        26	     ntotal=ntotal+1
                 if(ntotal .gt. nt-50) nsatura=1
                 npa(ntotal)=naux(j)
                 nmmm=nmmm+1
                 xcet(ntotal)=x(npp)
                 ycet(ntotal)=y(npp)
                 dxx(ntotal)=0.d0
                 particle(ntotal)=.true.
                 nnod(ntotal)=npointer(npp)
                 npart(npointer(npp))=ntotal

        25       continue
        28       continue
               end do
             enddo
             nh2(naux(j))=ntotal
             nh1(naux(j))=ntotal-nmmm+1
           end do
           if(nm .ne. 0) goto 6

           !Check that all particles have been found
           ncount=0
           do i=1,ntotal
             if(particle(i)) ncount=ncount+1
           enddo
           write(2,'(2x,6(2x,i9))') l,kl,n,ncount,ntotal,nsatura

           if(profiling) call profile(1,'   tree_calc')
           ! *** initial setting of h with ni0r.
           if(kl.eq.1) then
             do i=1,n11
               if(nv(i) .le. 20) then
                 nv(i)=20
!                 hfac=15.d0
!                 hexp=1.d0/4.d0
                 hfac=31.d0
                 hexp=1.d0/5.d0
               else if(nv(i) .ge. 150) then
                 nv(i)=150
                 hfac=31.d0
                 hexp=1.d0/5.d0
               else
                 hfac=31.d0
                 hexp=1.d0/5.d0
               endif
!                hfac=15.d0
!                hexp=1.d0/4.d0
               h0=h(i)
               h(i)=h(i)*0.5d0*(1.d0+hfac*dble(ni0r)/dble(nv(i)))**hexp
               h(i)=0.5d0*(h0+h(i))
             enddo
           endif

           if(profiling) call profile(0,'')
           if(ncount.eq.n .or. l.le.100) then
               !!BUILDING THE LISTED-LINK of NEIGHBORS
               !!Real and axis-reflective ghosts particles
               nv(1:n11)=0

               !$omp parallel private(i,nvSum,nodo,d1,d2,dc,j,d05,v1,v2)
               !$omp do schedule(guided)
               do i = 1,n11
 412              nvSum=0
                  nodo=0
411               nodo=nodo+1
                  d1=abs(a(i,1)-xcet(nodo))
                  d2=abs(a(i,2)-ycet(nodo))
                  dc=4.d0*h(i)+dxx(nodo)/2.d0
                  !     testea la interseccion del cubo jerarquico con el cubo de
                  !     accion de la particula "i".

                  if((d1.lt.dc).and.(d2.lt.dc))goto 451

421               if(nodo.eq.nh2(npa(nodo))) then
                     nodo=npa(nodo)
                     if(nodo.eq.0) goto 4011
                     goto 421
                  else
                     goto 411
                  endif
                  ! 451               if(dxx(nodo).eq.0.d0) then
451               if(particle(nodo)) then
                     j=nnod(nodo)
                     if(i.eq.j) goto 421
                     d05=sqrt(d1**2+d2**2)
                     v1=d05/h(i)
                     v2=d05/h(j)
                     if((v1.ne.0.d0.and.v1.lt.2.d0).or.(v2.lt.2.d0)) then
!                        if(v1.ne.0.d0.and.v1.lt.2.d0) then
                       nvSum=nvSum+1
                       if(nvSum.le.(vecMax-20)) then
                         vecMatrix(i,nvSum)=j
                       else
!!!                          print *,"Too many neighbors",nvSum,i
!                          stop
                       endif
                     endif
                     goto 421
                   endif
                   nodo=nh1(nodo)-1
                   goto 411
 4011              continue
                   nv(i)=min(nvSum,vecMax)
                   !if number of neighbors lesser than 10 change hi
                   if(nv(i) .le. 10) then
                     nv(i)=0
                     h(i)=1.2d0*h(i)
                     goto 412
                   endif
                   !if number of neighbors is larger than 120 change hi
                   if(nv(i) .ge. 120) then
                     nv(i)=0
                     h(i)=0.9d0*h(i)
                     goto 412
                   endif
                end do
                !$omp end do
                !$omp end parallel
              endif      ! close if(ncount.eq.n)
              if(profiling) call profile(1,'  find_neigh')

              if(profiling) call profile(0,'')
              nvmax=maxval(nv)
              nvmin=minval(nv)
              invmax=maxloc(nv,DIM=1)
              invmin=minloc(nv,DIM=1)
              ro3dmax=maxval(ro3d)
              ro3dmin=minval(ro3d)
              irhomax=maxloc(ro3d,DIM=1)
              irhomin=minloc(ro3d,DIM=1)
              !Density calculation
              ro2d=0.d0
              !$omp parallel private(i,h1,jVec,j,d1,d2,d02,d05,v1,pi2v1,w1d,dter,signature)
              !$omp do
              do i = 1,n11
                h1=1.d0/h(i)**2
                do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  d1=a(i,1)-a(j,1)
                  d2=a(i,2)-a(j,2)
                  d02=d1**2+d2**2
                  d05=sqrt(d02)
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  if(v1.lt.2.d0) then
                    call func(pi2v1,w1d)
                  else
                    w1d=0.d0
                  endif
                  w1d=w1d**hindex
                  dter=h1*w1d
                  signature=1
                  if(axis(j)) signature=-1
                  ro2d(i)=ro2d(i)+mass(j)*dter*signature
                enddo
                !Self-contribution and VEs
                ro2d(i)=pkv*ro2d(i)+pkv*mass(i)*h1
                ro3d(i)=ro2d(i)/(2.d0*pi*(a(i,1)-desplr))
                vol(i)=mass(i)/ro2d(i)
              end do
              !$omp end do
              !$omp end parallel

              if(profiling) call profile(1,'density_calc')
              !EOS call (id. gas)
              do i=1,n11
               p(i)=rgas*ro3d(i)*temp(i)/mui+arad/3.d0*temp(i)**4
!               u(i)=1.5*rgas*temp(i)/mui+arad*temp(i)**4/ro3d(i)
               c(i)=dsqrt(rgas*temp(i)/mui)
              enddo
              !Alfven velocity
              do i=1,n11
                valfven2(i)=(B(i,1)**2+B(i,2)**2+B(i,3)**2)/mu0/ro3d(i)
                va_sroot(i)=dsqrt(c(i)**2+valfven2(i))
              enddo

              !calculate the Stress tensor Sij  (including Magnetic field)

              if(profiling) call profile(0,'')
              !$omp parallel private(k,Bmod2,i,j)
              !$omp do
              do k=1,n11
                Bmod2=(B(k,1)**2+B(k,2)**2+B(k,3)**2)
                B2mu0(k)=Bmod2/(2.d0*mu0)
                pro(k)=(p(k)+B2mu0(k))/ro2d(k)
                proener(k)=p(k)/ro2d(k)
                do i=1,3
                  do j=1,3
                    Sij(k,i,j)=B(k,i)*B(k,j)/mu0/ro2d(k)
                  enddo
                enddo
              enddo
              !$omp end do

              !re-assign magnitudes (real+ ghosts).

              !$omp do
              do i=1,n
                h(i)=h(ngh(i))
                ro2d(i)=ro2d(ngh(i))
                ro3d(i)=ro3d(ngh(i))
                pro(i)=pro(ngh(i))
                c(i)=c(ngh(i))
                u(i)=u(ngh(i))
                p(i)=p(ngh(i))
                vol(i)=vol(ngh(i))
                valfven2(i)=valfven2(ngh(i))
                va_sroot(i)=va_sroot(ngh(i))
                do j=1,3
                  do k=1,3
                    if(i.eq.2 .and. j.eq.3 .and. axis(j)) then
                      Sij(i,j,k)=-Sij(ngh(i),j,k)
                    else
                      Sij(i,j,k)=Sij(ngh(i),j,k)
                    endif
                  enddo
                enddo
              enddo
              !$omp end do
              !$omp end parallel

              if(profiling) call profile(1,' Stress_tens')

              !IAD matrix calculation
              if(profiling) call profile(0,'')
              checkInorm=0.d0
              checkdeltax=0.d0
              checkdeltay=0.d0

              !$omp parallel private(i,hhi,h1,txx,txy,tyy,jVec,j,d1,d2,d02,d05,&
              !$omp    v1,pi2v1,w1d,dter,val)
              !$omp do
              do i = 1,n11
                hhi=h(i)*h(i)
                h1=1.d0/hhi
                txx=0.d0
                tyy=0.d0
                txy=0.d0
                do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  d1=abs(a(i,1)-a(j,1))
                  d2=abs(a(i,2)-a(j,2))
                  d02=d1**2+d2**2
                  d05=sqrt(d02)
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  if(v1.lt.2.d0) then
                    call func(pi2v1,w1d)
                  else
                    w1d=0.d0
                  endif
                  w1d=pkv*w1d**hindex

                  !calculation of  matrix Tau
                  dter=h1*w1d
                  val=vol(j)*dter
                  txx=txx+val*(a(j,1)-a(i,1))**2
                  tyy=tyy+val*(a(j,2)-a(i,2))**2
                  txy=txy+val*(a(j,1)-a(i,1))*(a(j,2)-a(i,2))
                  checkInorm(i)=checkInorm(i)+val
                  checkdeltax(i)=checkdeltax(i)+ (a(j,1)-a(i,1))*val
                  checkdeltay(i)= checkdeltay(i)+(a(j,2)-a(i,2))*val
                enddo

                !Matrix inversion
                cxx(i)=txx-txy**2/tyy
                cxx(i)=1.d0/cxx(i)
                cxy(i)=-cxx(i)*txy/tyy
                cyy(i)=tyy-txy**2/txx
                cyy(i)=1.d0/cyy(i)

                !Self-contribution to checkInorm
                checkInorm(i)=checkInorm(i)+mass(i)/ro2d(i)/h(i)**2*pkv
              end do
              !$omp end do
              !$omp end parallel

              !Matrix elements (including ghosts)
              do i=1,n
                cxx(i)=cxx(ngh(i))
                cxy(i)=cxy(ngh(i))
                cyy(i)=cyy(ngh(i))
              enddo
          if(profiling) call profile(0,'')
!         Calculation of gravity
           !$omp parallel private(i,rr2,rr,rr2_32)
           !$omp do
         do i=1,n11
           rr2=(a(i,1)-desplr)**2+(a(i,2)-desplz)**2
           rr=sqrt(rr2)
           ux=(a(i,1)-desplr)/rr
           uy=(a(i,2)-desplz)/rr
!           rr2_32=(rr2+rs**2)**1.5
           rr2=(rr2+rs**2)
           fgravr(i)=-g*masaPSR/rr2*ux
           rr2=rr**2
           fgravz(i)=-g*masaPSR/rr2*uy
         enddo
           !$omp end do
           !$omp end parallel
         if(profiling) call profile(1,'   grav_calc')
!!!            if(l.eq.(initial_adjust_iterations+1) .and. kl.eq.1) then 
!!!               do i=1,n11
!                this is a try to balance the Z-movement tuning temp
!!!                dmy=-(fgravz(i)+sum2(i))*temp(i)/sum2(i) !dmy is correction in temp
!!!                dmy1=min(temp(i)*0.2,abs(dmy))
!!!                temp(i)=temp(i)+dmy1*(dmy/abs(dmy))
!!!                temp0(i)=temp(i)
!!!             u(i)=1.5*rgas*temp(i)/mui+arad*temp(i)**4/ro3d(i)
!!!             u0(i)=u(i)
!            ---------------------------------------------------
!!!               enddo 
!!!             endif
            if(l.eq.(initial_adjust_iterations+1) .and. kl.eq.1) then 
               do i=1,n11
               ! Warning. If program crashes at the very beginning then cancel next line
                dmy=fgravr(i)+sum1(i)
                if(dmy.lt.0.d0) then 
!                v(i,3)=sqrt(-(fgravr(i)+sum1(i))*(a(i,1)-desplr))
                 v(i,3)=sqrt(-dmy*(a(i,1)-desplr))
                else
                 v(i,3)=sqrt(+dmy*(a(i,1)-desplr))
                endif
                v0(i,3)=v(i,3)
                lz0(i)=v(i,3)*(a(i,1)-desplr)
                lz(i)=lz0(i)
                v(i,1)=0.d0
                v(i,2)=0.d0
                v0(i,1)=v(i,1)
                v0(i,2)=v(i,2)                
!            ---------------------------------------------------
               enddo 
             endif


              !     Summations to find the viscous stress tensor and divergence and curl
              !     of the velocity

              if(profiling) call profile(0,'')
              divv=0.d0
              curlz=0.d0
              divB=0.d0
              divvrr=0.d0
              divvrz=0.d0
              divvzr=0.d0
              divvzz=0.d0
              divvphir=0.d0
              divvphiz=0.d0
              gradrBr=0.d0
              gradrBphi=0.d0
              gradrBz=0.d0
              f=0.d0

              !$omp parallel private(i,hhi,h1,h11,jVec,j,djix,djiy,d1,d2,d02,d05,&
              !$omp   v1,pi2v1,w1d,dteri,vjix,vjiy,vjiz,kernxxi,kernxyi,kernyxi,kernyyi,&
              !$omp   gradxha,gradyha,signature,dvx,dvy,dmy,dmy1,dmy2,smdiv)
              !$omp do
              do i = 1,n11
                hhi=h(i)*h(i)
                h1=1.d0/hhi
                h11=h1/hhi
                do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  djix=a(j,1)-a(i,1)
                  djiy=a(j,2)-a(i,2)
                  d1=abs(djix)
                  d2=abs(djiy)
                  d02=d1**2+d2**2
                  d05=sqrt(d02)

                  !Divergence and curl calculation
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  w1d=0.d0
                  if(v1.lt.2.d0) then
                    call func(pi2v1,w1d)
                  else
                    w1d=0.d0
                  endif
                  w1d=pkv*w1d**hindex
                  dteri=h1*w1d
                  vjix=v(j,1)-v(i,1)
                  vjiy=v(j,2)-v(i,2)
                  vjiz=v(j,3)-v(i,3)
                  kernxxi=cxx(i)*djix
                  kernxyi=cxy(i)*djiy
                  kernyyi=cyy(i)*djiy
                  kernyxi=cxy(i)*djix
                  gradxha=(kernxxi+kernxyi)*dteri
                  gradyha=(kernyxi+kernyyi)*dteri
                  signature=1
                  if(axis(j)) signature=-1
                  dvx=vjix*gradxha
                  dvy=vjiy*gradyha
                  if(derivphi1) then   
                   dmy=vol(j)*signature !choose phi=1 to perform derivatives
                  else
                   dmy=vol(j)*ro2d(j)*signature  !choose phi=ro2d to perform derivatives
                  endif
                  divv(i)=divv(i)+(dvx+dvy)*dmy
                  dmy1=gradxha*dmy
                  dmy2=gradyha*dmy
                  curlz(i)=curlz(i)+(vjiy*gradxha-vjix*gradyha)*dmy
                  divvrr(i)=divvrr(i)+vjix*dmy1
                  divvrz(i)=divvrz(i)+vjix*dmy2
                  divvzr(i)=divvzr(i)+vjiy*dmy1
                  divvzz(i)=divvzz(i)+vjiy*dmy2
                  divvphir(i)=divvphir(i)+vjiz*dmy1
                  divvphiz(i)=divvphiz(i)+vjiz*dmy2
                  divB(i)=divB(i)+(B(j,1)-B(i,1))*dmy1+(B(j,2)-B(i,2))*dmy2
                  gradrBr(i)=gradrBr(i)+(B(j,1)-B(i,1))*dmy1
                  gradrBz(i)=gradrBz(i)+(B(j,2)-B(i,2))*dmy1
                  gradrBphi(i)=gradrBphi(i)+(B(j,3)-B(i,3))*dmy1
                enddo
                !Divide by ro2(i) only if derivphi1 is true
               if(derivphi1 .eq. .false.) then 
                divv(i)=divv(i)/ro2d(i)
                divvrr(i)=divvrr(i)/ro2d(i)
                divvrz(i)=divvrz(i)/ro2d(i)
                divvzr(i)=divvzr(i)/ro2d(i)
                divvzz(i)=divvzz(i)/ro2d(i)
                divvphir(i)=divvphir(i)/ro2d(i)
                divvphiz(i)=divvphiz(i)/ro2d(i)
                divB(i)=divB(i)/ro2d(i)
                curlz(i)=curlz(i)/ro2d(i)
                gradrBr(i)=gradrBr(i)/ro2d(i)
                gradrBz(i)=gradrBz(i)/ro2d(i)
                gradrBphi(i)=gradrBphi(i)/ro2d(i)
               endif
                !Add Hoop-Stress terms to divB,divv  ((smoothed close to the axis)
                smdiv=a(i,1)-desplr+epsdiv*h(i)
                divv(i)=divv(i)+v(i,1)/smdiv
                divB(i)=divB(i)+B(i,1)/smdiv
!!!!!!                divvphir(i)=-0.5/smdiv*v(i,3)  !A try with analytical derivative on dv_phy_dr
              end do
              !$omp end do
              !$omp end parallel

              !set values for ghosts
              do i=1,n
                divv(i)=divv(ngh(i))
                curlz(i)=curlz(ngh(i))
              enddo

              if(profiling) call profile(1,' Visc_tensor')
              !To clean or not to clean the divergence. That is the question ..

              if(cleandiv) then
                !Calculation of the derivative of parameter phidiv necessary to
                !clean the divergence of B with the recipe of Wissing et al (2020).
                if(profiling) call profile(0,'')

                !$omp parallel private(i,tau_a)
                !$omp do
                do i=1,n11
                  cha(i)=fclean*va_sroot(i)
                  if(l.eq.lstart) cha0(i)=cha(i)
                  tau_a=h(i)/cha0(i)/sigma_c
                  phidiv(i)=(phidiv0(i)/cha0(i))*cha(i)-(cha(i)*cha0(i)*divB(i)+&
                  &    phidiv0(i)/2.d0*divv(i)*cha(i)/cha0(i)+&
                  &    phidiv0(i)/tau_a*cha(i)/cha0(i))*dt
                enddo
                !$omp end do
                !$omp end parallel

                if(profiling) call profile(1,'   Clean_div')
              endif

              !set values for ghosts
              do i=1,n
                phidiv(i)=phidiv(ngh(i))
              enddo
!     *********************************************************
!     calculation of the different terms in the momentum equation
!     **********************************************************

              if(profiling) call profile(0,'')
              sum1=0.d0; sum2=0.d0; sum3=0.d0;
              sum1B=0.d0; sum2B=0.d0; sum3B=0.d0
              sum5=0.d0
              Bindis=0.d0; Bdivclean=0.d0
              maxvsignal=0.d0;
              evisco=0.d0; eviscoB=0.d0; eviscu=0.d0
              eviscoBbis=0.d0; Jdis_r=0.d0; Jdis_z=0.d0; Jdis_phi=0.d0

              !$omp parallel private(i,beta,hhi,h1,jVec,j,dmy1,dmy2,d1,d2,d02,hhj,dmyr,dmyz,dmyphy,&
              !$omp     h2,d05,uijx,uijy,v1,v2,pi2v1,pi2v2,w1d,w2d,djix,djiy,dteri,dterj,&
              !$omp     signatureC,signatureNC,kernxxi,kernxyi,kernyxi,kernyyi,kernxxj,&
              !$omp     kernxyj,kernyxj,kernyyj,gradxha,gradxhb,gradyha,gradyhb,di,dj,&
              !$omp     vij1,vij2,vij3,vijrij,rij,fij,wij,vijsignal,vijsignalu,ro3dij,&
              !$omp     divi,divj,roti,rotj,fbalsi,fbalsj,qiji,qijiu,qijj,qijju,proi,proj,&
              !$omp     Sij_i11,Sij_i12,Sij_i13,Sij_i21,Sij_i22,Sij_i23,Sij_i31,Sij_i32,Sij_i33,&
              !$omp     Sij_j11,Sij_j12,Sij_j13,Sij_j21,Sij_j22,Sij_j23,Sij_j31,Sij_j32,Sij_j33,&
              !$omp     proeneri,momentumxxp,momentumyyp,momentumxx,momentumyy,momentumzz,&
              !$omp     dmy01,dmy02,dmy11,dmy12,aviscox,aviscoy,dpxavu,dpyavu,&
              !$omp     Bij2,Fabha,Fabhb,vijsignalBi,vijsignalBj,vijsignaldis,viscBi,viscBj,&
              !$omp     viscBdis,dmy,dmyB,dmyBclean,tqij,smdiv,Bdist,gradxhab,gradyhab,vijsignalab,resistiv_ab,momentumxxd,momentumyyd,momentumzzd,&
              !$omp     momentumxxl,momentumyyl,momentumzzl,momentumxxa,momentumyya,momentumzza)
              !$omp do
              do i = 1,n11
                if(B2mu0(i).ne.0.d0) then
                  beta=p(i)/B2mu0(i)
                  phi(i)=0.d0
                  if(beta.lt.1.d0) phi(i)=2.d0
                  if(1.d0.le.beta .and. beta.le.2.d0) phi(i)=2.d0*(2.d0-beta)
                else
                  phi(i)=0.d0
                endif
                hhi=h(i)*h(i)
                h1=1.d0/hhi
                do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  dmy1=a(i,1)-a(j,1)
                  dmy2=a(i,2)-a(j,2)
                  d1=abs(dmy1)
                  d2=abs(dmy2)
                  d02=d1**2+d2**2
                  hhj=h(j)*h(j)
                  h2=1.d0/hhj
                  d05=sqrt(d02)
                  uijx=dmy1/d05
                  uijy=dmy2/d05
                  v1=d05/h(i)
                  v2=d05/h(j)
                  pi2v1=pi2*v1
                  pi2v2=pi2*v2
                  if(v1.lt.2.d0) then
                    call func(pi2v1,w1d)
                  else
                    w1d=0.d0
                  endif
                  if(v2.lt.2.d0) then
                    call func(pi2v2,w2d)
                  else
                    w2d=0.d0
                  endif

                  !(j)-(i) is ihat is used in IAD gradient formula.
                  djix=-dmy1
                  djiy=-dmy2
                  w1d=pkv*w1d**hindex
                  w2d=pkv*w2d**hindex
                  dteri=h1*w1d
                  dterj=h2*w2d
                  signatureC=1
                  signatureNC=1
                  if(axis(j)) then
                    if(crossedMom) then
                      signatureC=-1
                    else
                      signatureNC=-1
                    endif
                  endif
                  kernxxi=cxx(i)*djix*signatureNC
                  kernxyi=cxy(i)*djiy*signatureNC
                  kernxxj=cxx(j)*djix*signatureC
                  kernxyj=cxy(j)*djiy*signatureC
                  kernyyi=cyy(i)*djiy*signatureNC
                  kernyxi=cxy(i)*djix*signatureNC
                  kernyyj=cyy(j)*djiy*signatureC
                  kernyxj=cxy(j)*djix*signatureC
                  gradxha=(kernxxi+kernxyi)*dteri
                  gradxhb=(kernxxj+kernxyj)*dterj
                  gradyha=(kernyxi+kernyyi)*dteri
                  gradyhb=(kernyxj+kernyyj)*dterj
                  gradxhab=0.5*(gradxha+gradxhb)
                  gradyhab=0.5*(gradyha+gradyhb)

                  ! Artificial viscosity calculation. Includes Balsara limiters.

                  di=a(i,1)-desplr
                  dj=a(j,1)-desplr
                  vij1=v(i,1)-v(j,1)
                  vij2=v(i,2)-v(j,2)
                  vij3=v(i,3)-v(j,3)
                  vijrij=vij1*(a(i,1)-a(j,1))+vij2*(a(i,2)-a(j,2))
                  rij=sqrt((a(j,1)-a(i,1))**2+(a(j,2)-a(i,2))**2)
                  fij=1.d0
                  if(vijrij.lt.0.d0) then
                    wij=vijrij/rij
                    !vijsignal=alfaAV*(va_sroot(i)+va_sroot(j))/2.d0-betaAV*wij
                    vijsignal=alfaAV*(va_sroot(i)+va_sroot(j))/2.d0-betaAV*wij
                  ro3dij=0.5d0*(ro3d(i)+ro3d(j))
                    !vijsignalu=sqrt(abs(p(i)-p(j))/ro3dij)
                    vijsignalu=abs(wij)
                    if(vijsignal.lt.0.d0) vijsignal=0.d0
                    maxvsignal(i)=max(maxvsignal(i),vijsignal)
                    divi=abs(divv(i))
                    divj=abs(divv(j))
                    roti=abs(curlz(i))
                    rotj=abs(curlz(j))
                    if(balsara_lim .and.l.ge.initial_adjust_iterations) then
                      fbalsi=divi/(divi+roti+1.d-5*c(i)/h(i))
                      fbalsj=divj/(divj+rotj+1.d-5*c(j)/h(j))
                      fij=max(0.05d0,0.5*(fbalsi+fbalsj))
                      fbalsi=fij
                      fbalsj=fij
                    else
                      fbalsi=1.d0
                      fbalsj=1.d0
                    endif
                    qiji=-vijsignal*wij/2.d0/pi
                    qijiu=alfau*vijsignalu*(u(i)-u(j))
                    qijj=qiji
                    qijju=qijiu
                    qiji=qiji*fbalsi
                    qijj=qijj*fbalsj
                  else
                    qiji=0.d0
                    qijj=0.d0
                    qijiu=0.d0
                    qijju=0.d0
                  endif        !closing the AV algorithm

                  if(.not.crossedMom) then
                    proi=pro(i)*di/ro2d(i)
                    proj=pro(j)*dj/ro2d(j)
                    proeneri=proener(i)*di/ro2d(i)
                    Sij_i11=Sij(i,1,1)/ro2d(i)
                    Sij_i12=Sij(i,1,2)/ro2d(i)
                    Sij_i13=Sij(i,1,3)/ro2d(i)
                    Sij_i21=Sij(i,2,1)/ro2d(i)
                    Sij_i22=Sij(i,2,2)/ro2d(i)
                    Sij_i23=Sij(i,2,3)/ro2d(i)
                    Sij_i31=Sij(i,3,1)/ro2d(i)
                    Sij_i32=Sij(i,3,2)/ro2d(i)
                    Sij_i33=Sij(i,3,3)/ro2d(i)
                    Sij_j11=Sij(j,1,1)/ro2d(j)
                    Sij_j12=Sij(j,1,2)/ro2d(j)
                    Sij_j13=Sij(j,1,3)/ro2d(j)
                    Sij_j21=Sij(j,2,1)/ro2d(j)
                    Sij_j22=Sij(j,2,2)/ro2d(j)
                    Sij_j23=Sij(j,2,3)/ro2d(j)
                    Sij_j31=Sij(j,3,1)/ro2d(j)
                    Sij_j32=Sij(j,3,2)/ro2d(j)
                    Sij_j33=Sij(j,3,3)/ro2d(j)
                  else
                    proi=pro(i)*di/ro2d(j)
                    proj=pro(j)*dj/ro2d(i)
                    proeneri=proener(i)*di/ro2d(j)
                    Sij_i11=Sij(i,1,1)/ro2d(j)
                    Sij_i12=Sij(i,1,2)/ro2d(j)
                    Sij_i13=Sij(i,1,3)/ro2d(j)
                    Sij_i21=Sij(i,2,1)/ro2d(j)
                    Sij_i22=Sij(i,2,2)/ro2d(j)
                    Sij_i23=Sij(i,2,3)/ro2d(j)
                    Sij_i31=Sij(i,3,1)/ro2d(j)
                    Sij_i32=Sij(i,3,2)/ro2d(j)
                    Sij_i33=Sij(i,3,3)/ro2d(j)
                    Sij_j11=Sij(j,1,1)/ro2d(i)
                    Sij_j12=Sij(j,1,2)/ro2d(i)
                    Sij_j13=Sij(j,1,3)/ro2d(i)
                    Sij_j21=Sij(j,2,1)/ro2d(i)
                    Sij_j22=Sij(j,2,2)/ro2d(i)
                    Sij_j23=Sij(j,2,3)/ro2d(i)
                    Sij_j31=Sij(j,3,1)/ro2d(i)
                    Sij_j32=Sij(j,3,2)/ro2d(i)
                    Sij_j33=Sij(j,3,3)/ro2d(i)
                  endif
                  !          Separate diagonal and Off-diagonal parts of Sij
                  momentumxxp=proi*gradxha+proj*gradxhb
                  momentumyyp=proi*gradyha+proj*gradyhb
                  momentumxx=Sij_i11*gradxhab*di+Sij_j11*gradxhab*dj+Sij_i12*gradyhab*di+Sij_j12*gradyhab*dj

                  momentumxxa=Sij_i11*gradxhab*di+Sij_j11*gradxhab*dj+Sij_i12*gradyhab*di+Sij_j12*gradyhab*dj
                  momentumyya=Sij_i21*gradxhab*di+Sij_j21*gradxhab*dj+Sij_i22*gradyhab*di+Sij_j22*gradyhab*dj
                  momentumzza=Sij_i31*gradxhab*di+Sij_j31*gradxhab*dj+Sij_i32*gradyhab*di+Sij_j32*gradyhab*dj
                  momentumxxd=-Sij_i11*gradxhab*di+Sij_j11*gradxhab*dj-&
                          &   Sij_i12*gradyhab*di+Sij_j12*gradyhab*dj
                  momentumyyd=-Sij_i21*gradxhab*di+Sij_j21*gradxhab*dj-&
                          &   Sij_i22*gradyhab*di+Sij_j22*gradyhab*dj
                  momentumzzd=-Sij_i31*gradxhab*di+Sij_j31*gradxhab*dj-&
                          &   Sij_i32*gradyhab*di+Sij_j32*gradyhab*dj

                  momentumxxl=mix*momentumxxa+(1.d0-mix)*momentumxxd
                  momentumyyl=mix*momentumyya+(1.d0-mix)*momentumyyd
                  momentumzzl=mix*momentumzza+(1.d0-mix)*momentumzzd
                   if(mix.ne.-1.d0) then
                  momentumxx=momentumxxa+(momentumxxl-momentumxxa)/2.d0*phi(i)
                  momentumyy=momentumyya+(momentumyyl-momentumyya)/2.d0*phi(i)
                  momentumzz=momentumzza+(momentumzzl-momentumzza)/2.d0*phi(i)
                 else
                  momentumxx=Sij_j11*gradxhab*dj+&
                          & Sij_j12*gradyhab*dj
                  momentumyy=Sij_j21*gradxhab*dj+&
                          & Sij_j22*gradyhab*dj
                  momentumzz=Sij_j31*gradxhab*dj+&
                          & Sij_j32*gradyhab*dj
                 endif

                  dmy01=vol(i)*mass(j)/mass(i)*qiji
                  dmy02=vol(j)*qijj
                  aviscox=0.5*(dmy01*gradxha+dmy02*gradxhb)
                  aviscoy=0.5*(dmy01*gradyha+dmy02*gradyhb)

                  dmy11=vol(i)*mass(j)/mass(i)*qijiu
                  dmy12=vol(j)*qijju
                  dpxavu=0.5d0*(dmy11*gradxha+dmy12*gradxhb)
                  dpyavu=0.5d0*(dmy11*gradyha+dmy12*gradyhb)

                  sum1(i)= sum1(i) + mass(j)*momentumxxp+aviscox
                  sum2(i)= sum2(i) + mass(j)*momentumyyp+aviscoy

                  sum1B(i)=sum1B(i)+mass(j)*momentumxx
                  sum2B(i)=sum2B(i)+mass(j)*momentumyy
                  sum3B(i)=sum3B(i)+mass(j)*momentumzz

                  Bij2=(B(i,1)-B(j,1))**2+(B(i,2)-B(j,2))**2+(B(i,3)-B(j,3))**2
                  Fabha=(-djix*gradxha-djiy*gradyha)/rij
                  Fabhb=(-djix*gradxhb-djiy*gradyhb)/rij

                  !Dissipation is always acting on B. Obtain (dB/dt)_diss  (cartesian part called Bindis calculated as in  Wissing&Shen2020):
                  vijsignalBi=sqrt(valfven2(i)) !magnetosonic velocity
                  vijsignalBj=sqrt(valfven2(j))
                  !*******************************************************
                  dmyr=-vij3*dmy2
                  dmyz=vij3*dmy1
                  dmyphy=vij1*dmy2-vij2*dmy1
!!!!                  vijsignaldis=2.d0*sqrt(dmyr**2+dmyz**2+dmyphy**2)/rij  !!A variant of the magnetic signal velocity in dissipation
                  !*******************************************************
                  vijsignaldis=(vijsignalBi+vijsignalBj)  !! signal velocity = alfven velocity
                  viscBdis=0.5d0*alfaB*vijsignaldis
!!                  viscBi=0.5*alfaB*vijsignaldis
!!                  viscBi=0.5d0*alfaB*vijsignalBi
!!                  viscBj=0.5d0*alfaB*vijsignalBj
                  dmy=vol(j)*viscBdis*(Fabha+Fabhb)/2.d0
                  dmyB=dmy
                  dmyBclean=vol(j)*(phidiv(i)+phidiv(j))
                  Bindis(i,1)=Bindis(i,1)+(B(i,1)-B(j,1))*dmyB
                  Bindis(i,2)=Bindis(i,2)+(B(i,2)-B(j,2))*dmyB
                  Bindis(i,3)=Bindis(i,3)+(B(i,3)-B(j,3))*dmyB
                  Bdivclean(i,1)=Bdivclean(i,1)+dmyBclean*gradxhab
                  Bdivclean(i,2)=Bdivclean(i,2)+dmyBclean*gradyhab
                  Bdivclean(i,3)=0.d0

                  dmyB=-dmy/ro3d(i)/mu0*Bij2  ! Apparently, there is a 1/2 constant missing here but it is placed in the energy Eq.

                  eviscoB(i)=eviscoB(i)+dmyB

                  tqij=mass(j)*(vij1*gradxha+vij2*gradyha)
                  sum5(i)=sum5(i)+tqij*proeneri
                  evisco(i)=evisco(i)+ aviscox*vij1+aviscoy*vij2
                  eviscu(i)=eviscu(i)+dpxavu*uijx+dpyavu*uijy
! Alternative  algorithm to implement the contribution to energy of dissipation
! (using currents J_dis)

                vijsignalab=vijsignaldis/2.d0
           resistiv_ab=0.5*alfaB*vijsignalab*rij
           resistiv_ab=sqrt(resistiv_ab)
            Jdis_r(i)=Jdis_r(i)+vol(j)*resistiv_ab*&
           &   (B(j,3)-B(i,3))*gradyhab
             Jdis_z(i)=Jdis_z(i)+vol(j)*resistiv_ab*&
           &   (B(j,3)-B(i,3))*gradxhab
             Jdis_phi(i)=Jdis_phi(i)+vol(j)*resistiv_ab*&
           &   ((B(j,1)-B(i,1))*gradyhab-&
           &   ((B(j,2)-B(i,2))*gradxhab))
                enddo   ! Close momentum neighbors loop

                !add non-axial part to Bindis
                smdiv=epsdiv*h(i)
                Bdist=min(di,h(i))
                viscBi=0.5d0*alfaB*vijsignalBi
                Bindis(i,1)=Bindis(i,1)+viscBi*Bdist*(gradrBr(i)-B(i,1)/di)/(di+smdiv)
                Bindis(i,2)=Bindis(i,2)+viscBi*Bdist*(gradrBz(i)/(di+smdiv))
                Bindis(i,3)=Bindis(i,3)+viscBi*Bdist*(gradrBphi(i)-B(i,3)/di)/(di+smdiv)
                Jdis_z(i)=Jdis_z(i)+sqrt(viscBi*Bdist)*B(i,3)/(di+smdiv)
                Jdis_r(i)=-Jdis_r(i)
              enddo !Close main momentum loop
              !$omp end do
              !$omp end parallel

              if(profiling) call profile(1,'    mom_ener')

               do i=1,n11
                  Jdis_r(i)=Jdis_r(i)/mu0
                  Jdis_z(i)=Jdis_z(i)/mu0
                  Jdis_phi(i)=Jdis_phi(i)/mu0
                  eviscoBbis(i)=(Jdis_r(i)**2+Jdis_z(i)**2+&
               & Jdis_phi(i)**2)
                   eviscoBbis(i)=eviscoBbis(i)/ro3d(i)
                enddo


              !         HOOP-STRESS TERMS
              if(profiling) call profile(0,'')
              !$omp parallel private(i,di,Btot2,hoop_r,hoop_phi,hoop)
              !$omp do
              do i=1,n11
                di=a(i,1)-desplr
                !with magnetic field in the Hoop-stress
                Btot2=B(i,1)**2+B(i,2)**2+B(i,3)**2
                hoop_r=2.d0*pi*((p(i)+Btot2/2.d0/mu0)-B(i,3)**2/mu0)/ro2d(i)
                sum1(i)=hoop_r-2.d0*pi*sum1(i)+2.d0*pi*sum1B(i)
                sum2(i)=-2.d0*pi*sum2(i)+2.d0*pi*sum2B(i)
                hoop_phi=2.d0*pi*(B(i,1)*B(i,3)/mu0)/ro2d(i)
                sum3(i)=2.d0*pi*sum3B(i)+hoop_phi
                evisco(i)=max(0.d0,evisco(i))
                hoop=-2.d0*2.d0*pi*p(i)/ro2d(i)*v(i,1)
                sum5(i)=hoop+2.d0*2.d0*pi*sum5(i)+ 2.d0*pi*evisco(i)+ 2.d0*eviscu(i)
!!               without AV contribution to energy equation
!!!               sum5(i)=hoop+2.d0*2.d0*pi*sum5(i) + 2.d0*eviscu(i)
              enddo
              !$omp end do
              !$omp end parallel
              if(profiling) call profile(1,' Hoop_stress')


              do i=1,n
                if(axis(i)) then
                  sum1(i)=-sum1(ngh(i))
                else
                  sum1(i)=sum1(ngh(i))
                endif
              enddo


              !           INDUCTION EQUATION
              if(profiling) call profile(0,'')
              !$omp parallel private(i,di,ratio,smdiv,Bcx,Bcz,Bcphi,lcentered,smdi,&
              !$omp    a11,a12,a13,a21,a22,a23,a31,a32,a33,dBr,dBtheta,dBphi,xdisk,dist,polar_ang)
              !$omp do
              do i=1,n11
               dist=sqrt((a(i,1)-desplr)**2+(a(i,2)-desplz)**2)
               polar_ang=atan(abs((a(i,1)-desplr)/(a(i,2)-desplz)))
               if((a(i,2)-desplz).lt.0.d0) polar_ang=pi-polar_ang
                xdisk=dist/rLC
!!!!               call derBdeutsch(psi,xdisk,polar_ang,rLC,Bp,P_PSR,rns,&
!!!!             &   omegaPSR,phi0,tt,dBr,dBtheta,dBphi)
                di=abs(a(i,1)-desplr)
                ratio=di/h(i)
                smdiv=epsdiv*h(i)
                a11=-(divvzz(i)+v(i,1)/(di+smdiv))
                a12=divvrz(i)
                a13=-v(i,3)/(di+smdiv)
                a21=divvzr(i)
                a22=-(divvrr(i)+v(i,1)/(di+smdiv))
                a23=0.d0
                a31=divvphir(i)
                a32=divvphiz(i)
                a33=-(divvrr(i)+divvzz(i))
                Bcx=B(i,1)
                Bcz=B(i,2)
                Bcphi=B(i,3)
                B(i,1)=B0(i,1)+(a11*Bcx+a12*Bcz+a13*Bcphi+ &
                 & Bindis(i,1))*dt-Bdivclean(i,1)*dt
                B(i,2)=B0(i,2)+(a21*Bcx+a22*Bcz+a23*Bcphi+ &
                 & Bindis(i,2))*dt-Bdivclean(i,2)*dt
                B(i,3)=B0(i,3)+(a31*Bcx+a32*Bcz+a33*Bcphi+&
                 & Bindis(i,3))*dt-Bdivclean(i,3)*dt
                lcentered=lz(i)
                smdi=di+epsdiv*h(i)
                v(i,1)=v0(i,1)+sum1(i)*dt+fgravr(i)*dt+lcentered**2/smdi**3*dt
                v(i,2)=v0(i,2)+sum2(i)*dt+fgravz(i)*dt
                lz(i)=lz0(i)+smdi*sum3(i)*dt
                v(i,3)=lz(i)/smdi
                u(i)=max(1.d-6, u0(i)+0.5d0*sum5(i)*dt+0.5*eviscoB(i)*dt)
!!!               u(i)=max(1.d-6, u0(i)+0.5d0*sum5(i)*dt+eviscoBbis(i)*dt)
                if(i.eq.1 .and.kl.eq.2) then 
                  write(9,'(12(1x,1e14.7))') tt,a(1,1)-desplr,a(1,2)-desplz,&
                          &v(1,1),v(1,2),v(1,3),&
                          &B(1,1),B(1,2),B(1,3),temp(1),ro3d(1),dble(nv(48))
                endif
              enddo
              !$omp end do
              !$omp end parallel
              if(profiling) call profile(1,'   Induction')

!             CALCULATION OF TEMPERATURE NEW MODEL
               ready(:)=.false.
             if(profiling) call profile(0,'')
              !$omp parallel private(i,k,ueos,dudt1,gnr,dgnr,dtemp)
              !$omp do schedule(static)
            outer:do i=1,n11
               temp(i)=temp0(i)  !initial guess of temperature
             inner: do k=1,itmaxNR
             if(ready(i)) cycle outer
               ueos=1.5*rgas*temp(i)/mui+arad*temp(i)**4/ro3d(i)  
               dudt1=1.5*rgas/mui+4.d0*arad*temp(i)**3/ro3d(i)
               gnr=u(i)-ueos
               dgnr=-dudt1
               dtemp=-gnr/dgnr 
               temp(i)=temp(i)+dtemp
               if(abs(dtemp/temp(i)) .le. tolNR) ready(i)=.true.
             end do inner
           end do outer
            !$omp end do
              !$omp do
              do i=1,n11
               p(i)=rgas*ro3d(i)*temp(i)/mui+arad/3.d0*temp(i)**4
               c(i)=dsqrt(rgas*temp(i)/mui)
               enddo
            !omp end do
            !$omp end parallel
            if(profiling) call profile(1,' Newton-Raphson')
!            ***  END IMPLICIT SEARCH OF TEMP ***


              !make the average of radial velocity near the axis
              if(averVelr) then
                if(profiling) call profile(0,'')
                !$omp parallel private(i,di,ratio,jVec,j,d1,d2,d02,d05,v1)
                !$omp do
                do i=1,n11
                  di=a(i,1)-desplr
                  ratio=abs(di)/h(i)
                  nvvv(i)=1
                  avevel(i)=0.d0
                  bfieldr(i)=0.d0
                  if(ratio.le.2.d0) then
                    do jVec=1,nv(i)
                      j=vecMatrix(i,jVec)
                      d1=abs(a(i,1)-a(j,1))
                      d2=abs(a(i,2)-a(j,2))
                      d02=d1**2+d2**2
                      d05=sqrt(d02)
                      v1=d05/h(i)
                      if(v1.ne.0.d0.and.v1.lt.2.d0) then
                        nvvv(i)=nvvv(i)+1
                        avevel(i)=avevel(i)+v(j,1)
                        bfieldr(i)=bfieldr(i)+B(j,1)
                      endif
                    enddo
                    avevel(i)=avevel(i)+v(i,1)  !Self-contribution
                    bfieldr(i)=bfieldr(i)+B(i,1)
                   endif
                 enddo
                 !$omp end do

                 !$omp do
               do i=1,n11
                  di=a(i,1)-desplr
                  if(nvvv(i).gt.1) then
                   v(i,1)=avevel(i)/nvvv(i)
                   B(i,1)=bfieldr(i)/nvvv(i)
                   if(di/h(i).le.0.5) B(i,1)=0.d0
                  endif
                enddo
                !$omp end do
                !$omp end parallel
                if(profiling) call profile(1,' Aver_Bfield')

              endif
              do i=1,n11
                a(i,1) = a0(i,1) + ((v(i,1)+v0(i,1))/2.d0)*dt
                a(i,2) = a0(i,2) + ((v(i,2)+v0(i,2))/2.d0)*dt
                if(a(i,1)-desplr .lt.0.d0) then
                  print*, 'crossing the axis !!'
                  print*, l,kl,i
                  print*, a(i,1)-desplr,a(i,2)-desplz
                  print*, a0(i,1)-desplr,a0(i,2)-desplz,ro3d(i)
                  print*, u(i),p(i),c(i)
                  print*, u0(i),h(i)
                  print*, v(i,1),v(i,2)
                  print*, v0(i,1),v0(i,2)
                  print*, nv(i)
                  stop
                endif
              enddo


              !Apply PBC
!              do i=1,n11
!                if(a(i,2)-desplz .gt. rtzu) a(i,2)=a(i,2)-(rtzu+rtzd)
!                if(a(i,2)-desplz .lt. -rtzd) a(i,2)=a(i,2)+(rtzu+rtzd)
!              enddo

              do i=1,n22
                ngh(i)=i
                axis(i)=.false.
              enddo

              do jg=1,nc
                call ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,jg,n,&
                       &    a,v,B,n11,n22,axis,ngh)
                do i=1,n
                  h(i)=h(ngh(i))
                  mass(i)=mass(ngh(i))
                  u(i)=u(ngh(i))
                  ro2d(i)=ro2d(ngh(i))
                enddo
              enddo
        end do

        !     **************************************************************
        !      END OF MODEL CALCULATION
        !    *************************************************************
        do i=1,n11
          ro3d(i)=ro2d(i)/(2.d0*pi*(a(i,1)-desplr))
        enddo

        ! Center of mass ***

        ycm = 0.0d0
        dmass=0.d0
        do i = 1,n11
          dmass=dmass+mass(i)
          ycm = ycm + a(i,2)*mass(i)
        end do
        ycm_check=ycm/dmass
        ycm=desplz
        xcm=desplr

        !Find max density
        rhomax=maxval(ro3d, DIM=1)
        irhomax=maxloc(ro3d, DIM=1)

        !Energy conservation

        ekin=0.d0
        eint=0.d0
        emag=0.d0
        egrav=0.d0
        romax=-1.d0
        errdivBav=0.d0
        errdivBmax=0.d0
        errdivBmin=1.d20
        angmomr=0.d0
        angmomz=0.d0
        angmomphi=0.d0
        Bradav=0.d0
        Bthetaav=0.d0
        Bphiav=0.d0
        masstot=0.d0

        do i=1,n11
          angmomr=angmomr-mass(i)*lz0(i)/(a(i,1)-xcm)*(a(i,2)-ycm)
          angmomz=angmomz+mass(i)*lz0(i)
          angmomphi=angmomphi+mass(i)*((a(i,2)-ycm)*v(i,1)-(a(i,1)-xcm)*v(i,2))
!          vmod=(v(i,1)**2+v(i,2)**2+v(i,3)**2)
          vmod=v(i,1)**2+v(i,2)**2    !non-orbital velocity
          ekin=ekin+0.5d0*mass(i)*vmod
          eint=eint+mass(i)*u(i)
          Bmo12=dsqrt(B(i,1)**2+B(i,2)**2+B(i,3)**2)
          emag=emag+0.5d0*Bmo12**2/mu0*mass(i)/ro3d(i)
            if(Bmo12 .ne. 0.d0) then
              errdivB=h(i)*abs(divB(i))/Bmo12
              Bradav=Bradav+abs(B(i,1))*mass(i)
              Bthetaav=Bthetaav+abs(B(i,2))*mass(i)
              Bphiav=Bphiav+abs(B(i,3))*mass(i)
              masstot=masstot+mass(i)
            else
              errdivB=0.d0
            endif
            if(errdivB.ge.errdivBmax) errdivBmax=errdivB
            if(errdivB.le.errdivBmin) errdivBmin=errdivB
            errdivBav=errdivBav+errdivB*mass(i)
          if(ro3d(i).gt.romax) romax=ro3d(i)
        end do
             errdivBav=errdivBav/masstot
             Bradav=Bradav/masstot
             Bthetaav=Bthetaav/masstot
             Bphiav=Bphiav/masstot
        if(l.le.initial_adjust_iterations) then
          etot0=ekin+eint+emag
          angmomtot_0=dsqrt(angmomr**2+angmomz**2+angmomphi**2)
          cm0=abs(ycm_check-desplz)
        endif
        angmomtot=dsqrt(angmomr**2+angmomz**2+angmomphi**2)
        cm=abs(ycm_check-desplz)
        tempaver=(temp(part(1))+temp(part(2))+temp(part(3))+temp(part(4)) &
          & +temp(part(5)))/5.d0
         densaver=(ro3d(part(1))+ro3d(part(2))+ro3d(part(3))+ro3d(part(4)) &
          & +ro3d(part(5)))/5.d0
        tempaver2=(temp(part2(1))+temp(part2(2))+temp(part2(3))+temp(part2(4)) &
          & +temp(part2(5)))/5.d0
         densaver2=(ro3d(part2(1))+ro3d(part2(2))+ro3d(part2(3))+ro3d(part2(4)) &
          & +ro3d(part2(5)))/5.d0
        tempaver3=(temp(part3(1))+temp(part3(2))+temp(part3(3))+temp(part3(4)) &
          & +temp(part3(5)))/5.d0
         densaver3=(ro3d(part3(1))+ro3d(part3(2))+ro3d(part3(3))+ro3d(part3(4)) &
          & +ro3d(part3(5)))/5.d0
          if(l.le.initial_adjust_iterations) then
          etot0=ekin+eint+emag
          angmomtot_0=dsqrt(angmomr**2+angmomz**2+angmomphi**2)
          cm0=abs(ycm_check-desplz)
          tempaver_0=tempaver
          tempaver2_0=tempaver2
          tempaver3_0=tempaver3
          densaver_0=densaver
          densaver2_0=densaver2
          densaver3_0=densaver3
         endif
        open(4,file='PSR_DiskEnerx05B8T00_00',access='append')
        elost=(ekin+eint+emag)-etot0
        tempaver=(tempaver-tempaver_0)/tempaver_0
        tempaver2=(tempaver2-tempaver2_0)/tempaver2_0
        tempaver3=(tempaver3-tempaver3_0)/tempaver3_0
        densaver=(densaver-densaver_0)/densaver_0
        densaver2=(densaver2-densaver2_0)/densaver2_0
        densaver3=(densaver3-densaver3_0)/densaver3_0
        write(4,137)l,tt,ekin,eint,egrav,emag,elost/etot0,&
        &   errdivBmin,errdivBmax,errdivBav,(cm-cm0)/1.d8,&
        &  abs(angmomtot-angmomtot_0)/angmomtot_0,angmomtot_0,&
        &  angmomtot,romax,Bradav,Bthetaav,Bphiav,tempaver,densaver,&
        &tempaver2,densaver2,tempaver3,densaver3
        close(4)
137     format(2x,i8,30(2x,e16.9))


        !Timestep control
        ! In acceleration
        time_accel=1.d50
        do i = 1,n11
          accel=sqrt(sum1(i)**2+sum2(i)**2)
          if(accel.eq.0.d0) cycle
          dmy=C_accel*sqrt(h(i)/accel)
          time_accel=min(time_accel,dmy)
        end do

        ! In temperature
        time_temp=1.d50
        if(timeTemp) then
         do i = 1,n11
!          du=abs(u(i)-u0(i))
          du=abs(temp(i)-temp0(i))
          if(du.eq.0.d0) cycle
!          dmy=(u(i)*dtnma*tol_temp)/du
          dmy=(temp(i)*dtnma*tol_temp)/du
          time_temp=min(time_temp,dmy)
         end do
        endif

        ! Courant condition
        time_Cour=1.d50
        do i = 1,n11
          if(maxvsignal(i).gt.0.d0) then
            dmy=min(Cour*h(i)/va_sroot(i),Cour*h(i)/maxvsignal(i))
          else
            dmy=Cour*h(i)/va_sroot(i)
          endif
          time_Cour=min(time_Cour,dmy)
        end do

        ! In the axis-crossing time in convergent particles.
        time_axis=1.d50
        if(timeaxis) then
          do i = 1,n11
            if(v(i,1).ge.0.d0) cycle
            dmy=tol_taxis*(a(i,1)-desplr)/abs(v(i,1))
            time_axis=min(time_axis,dmy)
          end do
        endif

        ! Final time-step
        dtnme=dtnma
        dtnma = min(time_accel,time_temp,time_Cour,time_axis)
        if(dtnma .ge. 1.05d0*dtnme) dtnma=1.05d0*dtnme
        if(dtnma .le. 0.2d0*dtnme) dtnma=0.2d0*dtnme
!!!!        if(dtnma.lt.initial_timestep) dtnma=initial_timestep
        if(dtnma.ge.(1.*P_PSR/100.d0)) dtnma=1.*P_PSR/100.d0
        dtn=0.5d0*(dtnme+dtnma)
        tt=tt+dtnma

!-------------------------------------------------------------------
!       New sentences aimed at removing particles too close to the 
!       symmetry axis. The algorithm removes only one particle at each 
!       iteration and redefines the order of the remaining particles

          i_out=n11
          remove=.false.
        do i=1,n11 
         ax=a(i,1)-desplr
         ay=a(i,2)-desplz
         if(ax.le.xlim .or. (ay.gt.ylimpos .or. ay.lt.ylimneg)) then 
           remove=.true.
           i_out=i
           print*, 'particle removed',l,i
           print*, ax,ay
           exit
         endif 
        enddo 

      if(remove) then 
        do i=i_out+1,n11
          a(i-1,1)=a(i,1) 
          a(i-1,2)=a(i,2)
          B(i-1,1)=B(i,1)
          B(i-1,2)=B(i,2)
          B(i-1,3)=B(i,3)
          v(i-1,1)=v(i,1)
          v(i-1,2)=v(i,2)
          v(i-1,3)=v(i,3)
          u(i-1)=u(i)
          phidiv(i-1)=phidiv(i)
          cha(i-1)=cha(i)
          lz(i-1)=lz(i)
          temp(i-1)=temp(i)
          h(i-1)=h(i)
          nv(i-1)=nv(i)
         enddo
          n11=n11-1
          n22=n11*2
          do i=1,n22
             ngh(i)=i
             axis(i)=.false.
            enddo
          do jg=1,nc
            call ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,jg,n,&
               &    a,v,B,n11,n22,axis,ngh)
              do i=1,n
                h(i)=h(ngh(i))
                mass(i)=mass(ngh(i))
                u(i)=u(ngh(i))
               ro2d(i)=ro2d(ngh(i))
               enddo
             enddo
         endif   ! big loop form remove

! -------------------------------------------------------------------      
         



        ! Init for next time-step
        a0 = a
        v0 = v
        u0 = u
        B0 = B
        phidiv0 = phidiv
        cha0 = cha
        lz0 = lz
        temp0=temp

        ! *************************************************************
        !  Output
        ! *************************************************************

!         if(mod(l,iodelay).eq.0 .or. l.eq.20) then
         if(mod(l,iodelay).eq.0 .or. l.eq. 20) then
         prefix='pT00B800.'
         value1=idint(tt*1.e4)
!         value1=l
         write(suffix,'(i5.5)') value1
         nomfs1=prefix//suffix
       !call nomfils(nomfs1,nomfs2,value1,prefix)
         open(10,file=nomfs1)
         do i =1,n11
!	ax=a(i,1)-xcm
!        ay=a(i,2)-ycm
         ax=a(i,1)-desplr
         ay=a(i,2)-desplz
         rrr=sqrt(ax**2+ay**2)
         ux=ax/rrr
         uy=ay/rrr
         vr=v(i,1)*ux+v(i,2)*uy
         Jdis=sqrt(Jdis_r(i)**2+Jdis_z(i)**2+Jdis_phi(i)**2)
         dmy=p(i)*2.d0*mu0/(B(i,1)**2+B(i,2)**2+B(i,3)**2)
          write(10,'(32(1x,e19.12),2x,i4)')ax,ay,rrr,mass(i),& 
      &   h(i),c(i),temp(i),ro3d(i),&
      &   ro2d(i),v(i,1),v(i,2),v(i,3),u(i),p(i),&
      &   sum1(i),sum2(i),fgravr(i),fgravz(i),&
      &   B(i,1),B(i,2),B(i,3),dmy,divB(i),phi(i),lz(i),&
      & checkInorm(i),Jdis,Jdis_r(i),Jdis_z(i),Jdis_phi(i),&
      &     divvphir(i),divvphiz(i),nv(i)
          enddo
           close(10)
             write(6,*) '******',nomfs1
              write(6,*) tt
        endif

       if (nsatura .eq. 1)  then
	   open(11,file='excede')
	   write(11,*) nsatura
	   close(11)
!	   goto 
           close(3)
	endif

!    *************************************************************************

      enddo !main loop
      stop
    END PROGRAM


!          subroutina que ordena una lista de numeros mediante un
!          algoritmo Heapsort (Numerical Recipes 1989), indx(i) es la
!          posicion en la lista sin ordenar del numero que tiene
!          orden i en la nueva lista ordenada.



subroutine indexx (n,arrin,indx)
  !	   implicit real*8 (a-h,o-z)
  integer arrin,indx,q
  dimension arrin(n),indx(n)
  do 11 j=1,n
    indx(j)=j
    11	  continue
    if(n.eq.1) return
    l=n/2+1
    ir=n
    10  	    continue
    if(l.gt.1) then
      l=l-1
      indxt=indx(l)
      q=arrin(indxt)
    else
      indxt=indx(ir)
      q=arrin(indxt)
      indx(ir)=indx(1)
      ir=ir-1
      if(ir.eq.1) then
        indx(1)=indxt
        return
      endif
    endif
    i=l
    j=l+l
    20	     if(j.le.ir) then
    if(j.lt.ir) then
      if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1
    endif
    if (q.lt. arrin(indx(j))) then
      indx(i)=indx(j)
      i=j
      j=j+j
    else
      j=ir+1
    endif
    goto 20
  endif
  indx(i)=indxt
  goto 10
end

subroutine func(pi2x,w)
  IMPLICIT NONE
  double precision, intent(in)::pi2x
  double precision, intent(out)::w
  if(pi2x.eq.0.d0) then
    w=1.d0
  else
    w=sin(pi2x)/pi2x
  endif
  RETURN
END subroutine func


subroutine dfunc(pi2,pi2x,dw)
  IMPLICIT NONE
  double precision, intent(in):: pi2,pi2x
  double precision, intent(out)::dw
  if(pi2x.eq.0.d0)then
    dw=0.d0
  else
    dw=pi2*(1.d0/tan(pi2x)-1.d0/(pi2x))
  endif
  RETURN
END subroutine dfunc


!                                SUBROUTINE NOMFILS
!    --------------------------------------------------------------
!
      subroutine nomfils(nomfs,nomfs2,i,prefix)
      implicit none
      integer,INTENT(in)::i
      integer decmillar,millar,cent,dec,uni,rest
      character,INTENT(out)::nomfs*14,nomfs2*14
      character,INTENT(in)::prefix*6
      character sufix*5
      decmillar=i/10000+48
      rest=mod(i,10000)
      millar=rest/1000+48
      rest=mod(i,1000)
      cent=rest/100+48
      rest=mod(i,100)
      dec=rest/10+48
      uni=mod(i,10)+48
      sufix=char(decmillar)//char(millar)//char(cent)//char(dec)//&
       &    char(uni)
      nomfs=prefix//sufix
      nomfs2='sed.'//sufix
        return
      end

      subroutine ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,jg,ng,&
        &    a1,v1,B1,n1,nt,axis,ngh)
        implicit none

        integer i
        double precision dx,dy

        integer, intent(in)  :: jg,n1,nt
        integer, intent(out) :: ng
        integer, dimension(nt), intent(out) :: ngh

        double precision, intent(in) :: desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5

        logical, dimension(nt), intent(out) :: axis

        double precision, dimension(nt,2), intent(inout) :: a1
        double precision, dimension(nt,3), intent(inout) :: v1,B1

        !      box 1   axis-reflective ghosts
        if(jg.eq.1) then
          ng=n1
          do i=1,n1
            dx=a1(i,1)-desplr
            if(dx.le.rtr5) then
              ng=ng+1
              a1(ng,1)=-a1(i,1)+2.d0*desplr
              a1(ng,2)=a1(i,2)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=-v1(i,3)
              B1(ng,1)=-B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=-B1(i,3)
              axis(ng)=.true.
              ngh(ng)=i
            endif
          end do
        endif
        !           box 2  (periodic ghosts up)
        if(jg.eq.2) then
          do i=1,n1
            dy=a1(i,2)-desplz
            if(dy.ge.rtz5*4.d0) then
              ng=ng+1
              a1(ng,1)=a1(i,1)
              a1(ng,2)=a1(i,2)-(rtzu+rtzd)
              v1(ng,1)=v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=v1(i,3)
              B1(ng,1)=B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=B1(i,3)
              axis(ng)=.false.
              ngh(ng)=i
            endif
          end do
        endif
        !           box 3 (periodic ghosts down)
        if(jg.eq.3) then
          do i=1,n1
            dy=a1(i,2)-desplz
            if(dy.le.-rtz5*4.d0) then
              ng=ng+1
              a1(ng,1)=a1(i,1)
              a1(ng,2)=(rtzu+rtzd)+a1(i,2)
              v1(ng,1)=v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=v1(i,3)
              B1(ng,1)=B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=B1(i,3)
              axis(ng)=.false.
              ngh(ng)=i
            endif
          end do
        endif
        !           Box 4 (Reflective ghosts right.
        !                  Have to be reflective to preserve the mass of the particles)
        if(jg.eq.4) then
          do i=1,n1
            dx=a1(i,1)-desplr
            if(dx.ge.rtr5*4.d0) then
              ng=ng+1
              a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
              a1(ng,2)=a1(i,2)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=v1(i,3)
              B1(ng,1)=B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=B1(i,3)
              axis(ng)=.false.
              ngh(ng)=i
            endif
          end do
        endif
        !           Box 5 (periodic ghosts upper-left)
        if(jg.eq.5) then
          do i=1,n1
            dx=a1(i,1)-desplr
            dy=a1(i,2)-desplz
            if(dx.le.rtr5 .and. dy.ge.rtz5*4.d0) then
              ng=ng+1
              a1(ng,1)=-a1(i,1)+2.d0*desplr
              a1(ng,2)=a1(i,2)-(rtzu+rtzd)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=-v1(i,3)
              B1(ng,1)=-B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=-B1(i,3)
              axis(ng)=.true.
              ngh(ng)=i
            endif
          end do
        endif
        !           Box 6 (reflective ghosts upper-right. First periodic up then make reflective left.
        !                  (have to be reflective to preserve the mass of the particle))
        if(jg.eq.6) then
          do i=1,n1
            dx=a1(i,1)-desplr
            dy=a1(i,2)-desplz
            if(dx.ge.rtr5*4.d0 .and. dy.ge.rtz5*4.d0)then
              ng=ng+1
              a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
              a1(ng,2)=a1(i,2)-(rtzu+rtzd)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=v1(i,3)
              B1(ng,1)=B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=B1(i,3)
              axis(ng)=.false.
              ngh(ng)=i
            endif
          end do
        endif
        !        Box 7 (Periodic ghosts lower-left)
        if(jg.eq.7) then
          do i=1,n1
            dx=a1(i,1)-desplr
            dy=a1(i,2)-desplz
            if(dx.le.rtr5 .and. dy.le.-rtz5*4.d0) then
              ng=ng+1
              a1(ng,1)=-a1(i,1)+2.d0*desplr
              a1(ng,2)=(rtzu+rtzd)+a1(i,2)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=-v1(i,3)
              B1(ng,1)=-B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)=-B1(i,3)
              axis(ng)=.true.
              ngh(ng)=i
            endif
          end do
        endif
        !           Box 8 (Reflective ghosts lower-right. First periodic down. Then reflective left.)
        if(jg.eq.8) then
          do i=1,n1
            dx=a1(i,1)-desplr
            dy=a1(i,2)-desplz
            if(dx.ge.rtr5*4.d0 .and. dy.le.-rtz5*4.d0)then
              ng=ng+1
              a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
              a1(ng,2)=(rtzu+rtzd)+a1(i,2)
              v1(ng,1)=-v1(i,1)
              v1(ng,2)=v1(i,2)
              v1(ng,3)=v1(i,3)
              B1(ng,1)=B1(i,1)
              B1(ng,2)=B1(i,2)
              B1(ng,3)= B1(i,3)
              axis(ng)=.false.
              ngh(ng)=i
            endif
          end do
        endif
        return
      end subroutine

      subroutine gravity(n_grav_rows,tp,gi1,gi2,gi3)

        implicit none
        integer i
        integer, intent(in) :: n_grav_rows
        double precision, dimension(n_grav_rows), intent(out) :: tp,gi1,gi2,gi3

        open(33,file='I1I2.txt')
        do i=1,n_grav_rows
          read(33,*) tp(i),gi1(i),gi2(i),gi3(i)
        enddo
        close(33)
        return
      end subroutine gravity

     subroutine Bdeutsch(psi,xdisk,theta,R_cyl,B0,P_PSR,RNS,&
              &      omega,phi0,time,Brad,Btheta,Bphi)
      implicit none 
       double precision psi,xdisk,R_cyl,B0,P_PSR,omega,theta
       double precision alpha_B,RNS,alpha2,alpha3,alpha4,alpha5,alpha6
       double precision alpha2_1,denom,rr,rho,rho2,rho3,time
       double precision phi0,phi,phii,Brad,Btheta,Bphi
       double precision d1,d2,d3,d4,q1,q2,q3,q4
       double precision, parameter::pi=dacos(-1.d0)
       double precision, parameter::clight=3.0d10
           alpha_B=RNS*omega/clight
           alpha2=alpha_B**2
           alpha3=alpha_B**3
           alpha4=alpha_B**4
           alpha5=alpha_B**5
           alpha6=alpha_B**6
           alpha2_1=alpha_B**2+1.d0
           denom=alpha6-3.d0*alpha4+36.d0
           rr=xdisk*R_cyl
           rho=rr*omega/clight
           rho2=rho**2
           rho3=rho**3
           d1=(alpha_B*rho+1.d0)/alpha2_1
           d2=(rho-alpha_B)/alpha2_1
           d3=(1.d0+alpha_B*rho-rho2)/alpha2_1
           d4=((rho2-1.d0)*alpha_B+rho)/alpha2_1
           q1=(3.d0*rho*(6.d0*alpha3-alpha5)+(3.d0-rho2)*&
            &  (6.d0*alpha2-3.d0*alpha4))/denom
           q2=((3.d0-rho2)*(alpha5-6.d0*alpha3)+3.d0*rho*&
            &  (6.d0*alpha2-3.d0*alpha4))/denom
           q3=((rho3-6.d0*rho)*(alpha5-6.d0*alpha3)+&
            &  (6.d0-3.d0*rho2)*(6.d0*alpha2-3.d0*alpha4))/(rho2*denom)
           q4=((6.d0-3.d0*rho2)*(alpha5-6.d0*alpha3)+&
            &  (6.d0*rho-rho3)*(6.d0*alpha2-3.d0*alpha4))/(rho2*denom)
           phii=phi0-omega*time+rho-alpha_B
           Brad=2.d0*B0*RNS**3/rr**3*(dcos(psi)*dcos(theta)+&
                & dsin(psi)*dsin(theta)*(d1*dcos(phii)+d2*sin(phii)))
           Btheta=B0*RNS**3/rr**3*(dcos(psi)*dsin(theta)-&
         & dsin(psi)*dcos(theta)*((q1+d3)*dcos(phii)+(q2+d4)*sin(phii)))
           Bphi=B0*RNS**3/rr**3*dsin(psi)* &
               & (-(q2*dcos(2.d0*theta)+d4)*dcos(phii)+&
               & (q1*dcos(2.d0*theta)+d3)*dsin(phii))
           return
           end subroutine Bdeutsch


		
	subroutine derBdeutsch(psi,xdisk,theta,R_cyl,B0,P_PSR,RNS,&
		&      omega,phi0,time,dBr,dBtheta,dBphi)
	implicit none 
	double precision psi,xdisk,R_cyl,B0,P_PSR,omega,theta
        double precision alpha_B,RNS,alpha2,alpha3,alpha4,alpha5,alpha6
        double precision alpha2_1,denom,rr,rho,rho2,rho3,time
        double precision phi0,phi,phii,Brad,Btheta,Bphi,dBr,dBtheta,dBphi
        double precision d1,d2,d3,d4,q1,q2,q3,q4
        double precision, parameter::pi=dacos(-1.d0)
        double precision, parameter::clight=3.0d10
           alpha_B=RNS*omega/clight
           alpha2=alpha_B**2
           alpha3=alpha_B**3
           alpha4=alpha_B**4
           alpha5=alpha_B**5
           alpha6=alpha_B**6
           alpha2_1=alpha_B**2+1.d0
           denom=alpha6-3.d0*alpha4+36.d0
           rr=xdisk*R_cyl
           rho=rr*omega/clight
           rho2=rho**2
           rho3=rho**3
           d1=(alpha_B*rho+1.d0)/alpha2_1
           d2=(rho-alpha_B)/alpha2_1
           d3=(1.d0+alpha_B*rho-rho2)/alpha2_1
           d4=((rho2-1.d0)*alpha_B+rho)/alpha2_1
           q1=(3.d0*rho*(6.d0*alpha3-alpha5)+(3.d0-rho2)*&
            &  (6.d0*alpha2-3.d0*alpha4))/denom
           q2=((3.d0-rho2)*(alpha5-6.d0*alpha3)+3.d0*rho*&
            &  (6.d0*alpha2-3.d0*alpha4))/denom
           q3=((rho3-6.d0*rho)*(alpha5-6.d0*alpha3)+&
            &  (6.d0-3.d0*rho2)*(6.d0*alpha2-3.d0*alpha4))/(rho2*denom)
           q4=((6.d0-3.d0*rho2)*(alpha5-6.d0*alpha3)+&
            &  (6.d0*rho-rho3)*(6.d0*alpha2-3.d0*alpha4))/(rho2*denom)
           phi=phi0-omega*time
           phii=phi+rho-alpha_B
           dBr=2.d0*B0*RNS**3/rr**3*dsin(psi)*dsin(theta)*&
                & (-d1*dsin(phii)+d2*dcos(phii))*(-omega)
           dBtheta=-B0*RNS**3/rr**3*dsin(psi)*dcos(theta)*&
       &    (-(q1+d3)*dsin(phii)+(q2+d4)*dcos(phii))*(-omega)
           dBphi=B0*RNS**3/rr**3*dsin(psi)* &
               & ((q2*dcos(2.d0*theta)+d4)*dsin(phii)+&
               & (q1*dcos(2.d0*theta)+d3)*dcos(phii))*(-omega)
            dBr=0.d0
            dBtheta=0.d0
            dBphi=0.d0
           return
	end subroutine derBdeutsch


      subroutine profile(mode,label)
        integer, intent(in):: mode
        character :: label*12
        double precision start
        double precision omp_get_wtime
        double precision end

        save start

        if(mode.eq.0) then
          start=omp_get_wtime()
        else if(mode.eq.1) then
          end=omp_get_wtime()
          print '(a12,a2,es12.5,a2)',label(1:),': ',end-start,' s'
        endif
      end subroutine profile
