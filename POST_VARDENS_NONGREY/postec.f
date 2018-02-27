
      program postec

      use decomp_2d
      use decomp_2d_io
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'

      real mpi_cm(0:imax),mpi_sm(imax)
      real tmp1(0:imax),tmp2(0:imax),tmp3(0:imax),tmp4(0:imax),pow,Lz,rbulk(500),work(kmax)
      real utot(imax,jmax,kmax),vtot(imax,jmax,kmax),wtot(imax,jmax,kmax),ctot(imax,jmax,kmax)      
      real qptot(imax,jmax,kmax),Gptot(imax,jmax,kmax),kptot(imax,jmax,kmax)      
      real Gtot(imax,jmax,kmax),ktot(imax,jmax,kmax),vctot(imax,jmax,kmax)  
      real ptot(imax,jmax,kmax),qtot(imax,jmax,kmax),ytot(jmax),ztot(kmax),unew_s(imax,jmax/p_row,kmax/p_col)
      real utot1(imax,jmax,kmax),mpi_s(imax,jmax,kmax),ytot1(jmax),ztot1(kmax)
      real qstot(imax,jmax,kmax)
      real um(imax),vm(imax),wm(imax),pm(imax),cm(imax),qm(imax),Gm(imax),km(imax),snd(imax)
      real qpm(imax),Gpm(imax),kpm(imax),vcm(imax)
      real uft(imax,jmax,kmax),vft(imax,jmax,kmax)
      real qpft(imax,jmax,kmax),Gpft(imax,jmax,kmax)
      real kpft(imax,jmax,kmax)
      real vcft(imax,jmax,kmax)
      real wft(imax,jmax,kmax),cft(imax,jmax,kmax)
      real pft(imax,jmax,kmax),qft(imax,jmax,kmax)
      real Gft(imax,jmax,kmax),kft(imax,jmax,kmax)
      real hi_s(imax,4),hi_c(0:imax,4)
      real tmp5(0:imx,jmax,kmax/p_col)
      real tmp6(0:imx,jmax/p_col,kmax)
      real wspec(kmax)
      real spec(kmax/2)
       
      character*5 cha,cha2
      character*1 yes,ycor,yvort
      integer ire,nfiles,rnk,ier,loc
      integer dimens(2),coordr(2),coordn(2),coords(2),coordt(2),coordb(2)

      call mpi_init(ierr)
      call decomp_2d_init(imax+2,jmax,kmax,p_row,p_col)

      Rout=2.
      Rin=0.
      ytot=0
      ztot=0
      utot=0
      vtot=0
      wtot=0
      ctot=0
      ptot=0
      qtot=0
      qstot=0
      vctot=0
      Gtot=0
      ktot=0
      qptot=0
      Gptot=0
      kptot=0
      qpft=0
      Gpft=0
      kpft=0
      uft=0
      vft=0
      wft=0
      cft=0
      pft=0
      qft=0
      Gft=0
      kft=0
      rp=0
      open(5,file='INTEC')
      read(5,*)       iskip,Re,Lz,loc
            close(5)
      if (nrank.eq.0) then
      write(*,*)      iskip,Re,Lz
      endif
      
      re = int (RE)
      call cnvstr(ire,cha)
      call mkgrid(Lz,1,nrank)
      call loadd(0,ii,iskip)

      call inter_c_s_m(imax,jmax*kmax/px,unew,unew_s,dr)

        do k=1,kmax/p_col
         do j=1,jmax/p_row
           do i=1,imax
             um(i)=um(i)+unew_s(i,j,k)/(jmax*kmax)
             vm(i)=vm(i)+vnew(i,j,k)/(jmax*kmax)
             wm(i)=wm(i)+wnew(i,j,k)/(jmax*kmax)
             pm(i)=pm(i)+p(i,j,k)/(jmax*kmax)
             cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax)
             qm(i)=qm(i)+qrad(i,j,k)/(jmax*kmax)
             vcm(i)=vcm(i)+vcor(i,j,k)/(jmax*kmax)
             Gm(i)=Gm(i)+Gnew(i,j,k)/(jmax*kmax)
             km(i)=km(i)+knew(i,j,k)/(jmax*kmax)
             qpm(i)=qpm(i)+qpnew(i,j,k)/(jmax*kmax)
             Gpm(i)=Gpm(i)+Gpnew(i,j,k)/(jmax*kmax)
             kpm(i)=kpm(i)+kpnew(i,j,k)/(jmax*kmax)
            enddo
          enddo
        enddo
    
      call mpi_allreduce(um,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      um=snd
      call mpi_allreduce(vm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vm=snd
      call mpi_allreduce(wm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wm=snd
      call mpi_allreduce(cm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      cm=snd
      call mpi_allreduce(pm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      pm=snd
      call mpi_allreduce(qm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qm=snd
      call mpi_allreduce(vcm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vcm=snd
      call mpi_allreduce(Gm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gm=snd
      call mpi_allreduce(km,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      km=snd
      call mpi_allreduce(qpm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qpm=snd
      call mpi_allreduce(Gpm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gpm=snd
      call mpi_allreduce(kpm,snd,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kpm=snd


      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         jr=xstart(2)+j-1
         kr=xstart(3)+k-1
         ytot(jr)=(xstart(2)+j-2)*dtheta+dtheta/2.
         ztot(kr)=(xstart(3)+k-2)*dz+dz/2.                
         utot(i,jr,kr)=unew_s(i,j,k)
         vtot(i,jr,kr)=vnew(i,j,k)
         wtot(i,jr,kr)=wnew(i,j,k)
         ctot(i,jr,kr)=cnew(i,j,k)
         Gtot(i,jr,kr)=Gnew(i,j,k)
         ktot(i,jr,kr)=knew(i,j,k)
         ptot(i,jr,kr)=p(i,j,k)
         qtot(i,jr,kr)=qrad(i,j,k)
         qstot(i,jr,kr)=qsnew(i,j,k)
         vctot(i,jr,kr)=vcor(i,j,k)
         utot(i,jr,kr)=unew_s(i,j,k)
         qptot(i,jr,kr)=qpnew(i,j,k)
         Gptot(i,jr,kr)=Gpnew(i,j,k)
         kptot(i,jr,kr)=kpnew(i,j,k)
        end do
       end do
      end do
 
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         jr=xstart(2)+j-1
         kr=xstart(3)+k-1
         uft(i,jr,kr)=(unew_s(i,j,k)-um(i))
         vft(i,jr,kr)=(vnew(i,j,k)-vm(i))
         wft(i,jr,kr)=(wnew(i,j,k)-wm(i))
         cft(i,jr,kr)=(cnew(i,j,k)-cm(i))
         Gft(i,jr,kr)=(Gnew(i,j,k)-Gm(i))
         kft(i,jr,kr)=(knew(i,j,k)-km(i))
         pft(i,jr,kr)=(p(i,j,k)-pm(i))
         qft(i,jr,kr)=(qrad(i,j,k)-qm(i))
         vcft(i,jr,kr)=(vcor(i,j,k)-vcm(i))
         uft(i,jr,kr)=(unew_s(i,j,k)-um(i))
         qpft(i,jr,kr)=(qpnew(i,j,k)-qpm(i)) 
         Gpft(i,jr,kr)=(Gpnew(i,j,k)-Gpm(i))
         kpft(i,jr,kr)=(kpnew(i,j,k)-kpm(i))
        end do
       end do
      end do

      call mpi_allreduce(ytot,ytot1,jmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ytot=ytot1/p_col
      call mpi_allreduce(ztot,ztot1,kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ztot=ztot1/p_row
      call mpi_allreduce(utot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      utot=mpi_s
      call mpi_allreduce(vtot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vtot=mpi_s
      call mpi_allreduce(wtot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wtot=mpi_s
      call mpi_allreduce(ctot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ctot=mpi_s
      call mpi_allreduce(ptot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ptot=mpi_s
      call mpi_allreduce(qstot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qstot=mpi_s
      call mpi_allreduce(qtot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qtot=mpi_s
      call mpi_allreduce(vctot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vctot=mpi_s
      call mpi_allreduce(Gtot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gtot=mpi_s
      call mpi_allreduce(ktot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ktot=mpi_s
      call mpi_allreduce(qptot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qptot=mpi_s
      call mpi_allreduce(Gptot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gptot=mpi_s
      call mpi_allreduce(kptot,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kptot=mpi_s
      call mpi_allreduce(uft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      uft=mpi_s
      call mpi_allreduce(vft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vft=mpi_s
      call mpi_allreduce(wft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wft=mpi_s
      call mpi_allreduce(cft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      cft=mpi_s
      call mpi_allreduce(pft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      pft=mpi_s
      call mpi_allreduce(qft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qft=mpi_s
      call mpi_allreduce(vcft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vcft=mpi_s
      call mpi_allreduce(Gft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gft=mpi_s
      call mpi_allreduce(kft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kft=mpi_s
      call mpi_allreduce(qpft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qpft=mpi_s
      call mpi_allreduce(Gpft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gpft=mpi_s
      call mpi_allreduce(kpft,mpi_s,imax*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kpft=mpi_s

      if (nrank.eq.0) then
      open(45,file='tecplot/bench.dat')
      write(45,*) ' VARIABLES ="X","Y","Z", "U-vel","V-vel","W-vel","T","P","Q","Vc","G","K","Qp","Gp","kp","Qs"  '
      write(45,*) ' ZONE J= ',jmax, 'K = ',kmax ,'  F=POINT '
      do k=1,kmax
        do j=1,jmax
         write(45,'(20E16.6)') rp(loc),ytot(j),ztot(k),uft(loc,j,k),vft(loc,j,k),wft(loc,j,k),
     +               cft(loc,j,k),pft(loc,j,k),qft(loc,j,k),vcft(loc,j,k),Gft(loc,j,k),kft(loc,j,k)
     +               ,qpft(loc,j,k),Gpft(loc,j,k),kpft(loc,j,k),qstot(loc,j,k)
        end do
      end do
      close(45) 

      open(45,file='tecplot/benchtot.dat')
      write(45,*) ' VARIABLES ="X","Y","Z", "U-vel","V-vel","W-vel","T","P","Q","Vc","G","K","qp","Gp","kp","Qs"  '
      write(45,*) ' ZONE J= ',jmax, 'K = ',kmax ,'  F=POINT '
      do k=1,kmax
        do j=1,jmax
         write(45,'(20E16.6)') rp(loc),ytot(j),ztot(k),utot(loc,j,k),vtot(loc,j,k),wtot(loc,j,k),
     +               ctot(loc,j,k),ptot(loc,j,k),qtot(loc,j,k),vctot(loc,j,k),Gtot(loc,j,k),ktot(loc,j,k)
     +               ,qptot(loc,j,k),Gptot(loc,j,k),kptot(loc,j,k),qstot(loc,j,k)
        end do
      end do
      close(45) 

      open(45,file='tecplot/benchxzft.dat')
      write(45,*) ' VARIABLES ="X","Z", "U-vel","V-vel","W-vel","T","P","Q","Vc","G","K","qp","Gp","kp","Qs"  '
      write(45,*) ' ZONE I= ',imax, 'K = ',kmax ,'  F=POINT '
        do k=1,kmax
      do i=1,imax
         write(45,'(20E16.6)') rp(i),ztot(k),uft(i,jmax/2,k),vft(i,jmax/2,k),wft(i,jmax/2,k),
     +               cft(i,jmax/2,k),pft(i,jmax/2,k),qft(i,jmax/2,k),vcft(i,jmax/2,k),Gft(i,jmax/2,k),kft(i,jmax/2,k),
     +               qpft(i,jmax/2,k),Gpft(i,jmax/2,k),kpft(i,jmax/2,k),qstot(i,jmax/2,k)
        end do
      end do
      close(45) 

      open(45,file='tecplot/benchxztot.dat')
      write(45,*) ' VARIABLES ="X","Z", "U-vel","V-vel","W-vel","T","P","Q","Vc","G","K","qp","Gp","kp","Qs"  '
      write(45,*) ' ZONE I= ',imax, 'K = ',kmax ,'  F=POINT '
        do k=1,kmax
      do i=1,imax
         write(45,'(20E16.6)') rp(i),ztot(k),utot(i,jmax/2,k),vtot(i,jmax/2,k),wtot(i,jmax/2,k),
     +               ctot(i,jmax/2,k),ptot(i,jmax/2,k),qtot(i,jmax/2,k),vctot(i,jmax/2,k),Gtot(i,jmax/2,k),ktot(i,jmax/2,k),
     +               qptot(i,jmax/2,k),Gptot(i,jmax/2,k),kptot(i,jmax/2,k),qstot(i,jmax/2,k)
        end do
      end do
      close(45) 

!      open(45,file='tecplot/wholetot.dat')
!      write(45,*) ' VARIABLES ="X","Y","Z", "U-vel","V-vel","W-vel","T","P","Q","G","K"  '
!      write(45,*) ' ZONE I= ',imax, 'J = ',jmax,'K = ',kmax ,'  F=POINT '
!      do k=1,kmax
!      do j=1,jmax
!      do i=1,imax
!         write(45,'(20E16.6)') rp(i),ytot(j),ztot(k),utot(i,j,k),vtot(i,j,k),wtot(i,j,k),
!     +               ctot(i,j,k),ptot(i,j,k),qtot(i,j,k),Gtot(i,j,k),ktot(i,j,k)
!      end do
!      end do
!      end do
!      close(45) 
!
!      open(45,file='tecplot/wholeft.dat')
!      write(45,*) ' VARIABLES ="X","Y","Z", "U-vel","V-vel","W-vel","T","P","Q","G","K"  '
!      write(45,*) ' ZONE I= ',imax, 'J = ',jmax,'K = ',kmax ,'  F=POINT '
!      do k=1,kmax
!      do j=1,jmax
!      do i=1,imax
!         write(45,'(20E16.6)') rp(i),ytot(j),ztot(k),uft(i,j,k),vft(i,j,k),wft(i,j,k),
!     +               cft(i,j,k),pft(i,j,k),qft(i,j,k),Gft(i,j,k),kft(i,j,k)
!      end do
!      end do
!      end do
!      close(45) 
      endif

      call decomp_2d_finalize
      call mpi_finalize(ierr)
     
      end


      subroutine loadd(ini,nnrank,istap)
      use decomp_2d
      use decomp_2d_io
      include 'par_post.txt'
      include 'common.txt'
      real utmp(0:i1,jmax/p_row,kmax/p_col)
      real vtmp(0:i1,jmax/p_row,kmax/p_col)
      real wtmp(0:i1,jmax/p_row,kmax/p_col)
      real ptmp(0:i1,jmax/p_row,kmax/p_col)
      real ctmp(0:i1,jmax/p_row,kmax/p_col)

      integer istap
      character*5 cha
      character*5 cha2
      call cnvstr(nnrank,cha)
      call cnvstr(istap,cha2)
      ptmp = 0


      if (ini.eq.0) then
      if (nrank.eq.0)  write (*,*) cha2,imax,kmax,jmax
      call decomp_2d_read_one(1,unew,'DATAtex/u.'//cha2//'.dns')
      call decomp_2d_read_one(1,vnew,'DATAtex/v.'//cha2//'.dns')
      call decomp_2d_read_one(1,wnew,'DATAtex/w.'//cha2//'.dns')
      call decomp_2d_read_one(1,cnew,'DATAtex/c.'//cha2//'.dns')
      call decomp_2d_read_one(1,ptmp,'DATAtex/p.'//cha2//'.dns')
      call decomp_2d_read_one(1,qsnew,'DATAtex/qS.'//cha2//'.dns')
      call decomp_2d_read_one(1,qrad,'DATAtex/q.'//cha2//'.dns')
      call decomp_2d_read_one(1,vcor,'DATAtex/vc.'//cha2//'.dns')
      call decomp_2d_read_one(1,Gnew,'DATAtex/G.'//cha2//'.dns')
      call decomp_2d_read_one(1,knew,'DATAtex/k.'//cha2//'.dns')
      call decomp_2d_read_one(1,qpnew,'DATAtex/qp.'//cha2//'.dns')
      call decomp_2d_read_one(1,Gpnew,'DATAtex/Gp.'//cha2//'.dns')
      call decomp_2d_read_one(1,kpnew,'DATAtex/kp.'//cha2//'.dns')
      endif
        do i=1,imax
        p(i,:,:) = ptmp(i,:,:)
        enddo
       enew=cnew

      end



      subroutine mkgrid(Lz,uplus,rank)
      include 'par_post.txt'
      include 'common.txt'
      real const,Lz,rp_c(0:imax),rp_s(imax)
      dz =Lz/(kmax)
      dr = (Rout-Rin)/imax
      dtheta =6.*atan(1.)/jmax
       rmax = Rout-Rin
       ru(0)=0.

        do i=1,imax/2
          x  = 1.*i/imax
          dx = 0.5-1.45*(x-0.5)**2.!0.00145*x**3. -5.966*x**2. +734.37*x+1493.81
          ru(i)=ru(i-1)+dx
        enddo
        rnorm = ru(imax/2)
        do i=1,imax/2
          ru(i)=Rout/2.*ru(i)/rnorm
        enddo
       do i=imax,imax/2+1,-1
         ru(i)=Rout-ru(imax-i)
       enddo

      do i=1,imax
      rp(i)=0.5*(ru(i)+ru(i-1))
      enddo
      rp(0) = Rin - (rp(1)-Rin)
      rp(i1)=ru(imax)+(Ru(imax)-rp(imax))
      do i=1,imax
          rp_s(i)=rp(i)
      enddo
      call deriv_s_c(imax,rp_s,mr_c,dr)

      call inter_c_s(imax,mr_c,mr_s,dr)
      mr_c = 1./mr_c
      mr_s = 1./mr_s

      if (rank.eq.0) then
      open(11,file = 'grid.txt')
      write(11,*) Re,Ru(imax)
      do i=1,imax
         write(11,'(i5,4F12.6)') i,Ru(i),Rp(i)
      enddo
      endif
      close(11)


!      if (rank.eq.0) then
!      do i=1,imax
!      if (rank.eq.0)                write(6,111) i,(rp(i)-Rin),(Rout-rp(i)),(Ru(i)-Ru(i-1))
!      enddo
!      endif
!111   format ('Grid node =  ',i5, ' y+an =', F15.5, ' y+Pipe  =  ', F15.5, ' dy+  =  ', F15.5)
      end




      subroutine cnvstr(ing,cha)
      integer intg,i1,i2,i3,ing
      character*5 cha
      character c1,c2,c3,c4,c5
      intg=ing
      i1=intg/10000
      call cn(i1,c1)
      intg=intg-i1*10000
      i1=intg/1000
      call cn(i1,c2)
      intg=intg-i1*1000
      i1=intg/100
      call cn(i1,c3)
      intg=intg-i1*100
      i1=intg/10
      call cn(i1,c4)
      intg=intg-i1*10
      i1=intg
      call cn(i1,c5)
      cha=c1//c2//c3//c4//c5
      end

      subroutine cn(inte,ch)
      integer inte
      character ch
      if (inte .eq. 0) ch='0'
      if (inte .eq. 1) ch='1'
      if (inte .eq. 2) ch='2'
      if (inte .eq. 3) ch='3'
      if (inte. eq. 4) ch='4'
      if (inte. eq. 5) ch='5'
      if (inte. eq. 6) ch='6'
      if (inte. eq. 7) ch='7'
      if (inte. eq. 8) ch='8'
      if (inte. eq. 9) ch='9'
      end


