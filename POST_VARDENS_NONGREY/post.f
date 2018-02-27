
      program post

      use decomp_2d
      use decomp_2d_io
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'

      real um(0:imax),vm(imax),wm(imax),pm(imax),prr(imax),kkk(imax),km(imax),krr(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),dw(0:imax),dww(imax),w3(imax),w4(imax)
      real u3(imax),u4(imax),v3(imax),v4(imax),bulk1,bulkrms,ccc(imax),cbulk,www(imax)
      real mpi_cm(0:imax),mpi_sm(imax),opthk,stime2
      real uu(imax),vv(imax),ww(imax),duc_c(0:imax),duc(imax)
      real hic_s(imax,jmax/p_row,kmax/p_col),hic_c(0:imax,jmax/p_row,kmax/p_col)
      real tmp1(0:imax),tmp2(0:imax),tmp3(0:imax),tmp4(0:imax),pow,Lz,rbulk(500),work(kmax)
      real uMean(0:imax+1,jmax/p_row, kmax/p_col)
      real vMean(0:imax+1,jmax/p_row, kmax/p_col)
      real wMean(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc_tmp(0:imax+1,jmax/p_row, kmax/p_col)
      real vfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real rfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real pfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ducfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dvcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dwcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real cfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real qfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real kfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Gfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real efluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pefluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pafluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fxfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fyfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fzfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real umf(0:imax),vmf(imax),wmf(imax),vcm(imax)
      real urf(imax),vrf(imax),wrf(imax),uwf(imax)
      real uuf(imax),vvf(imax),wwf(imax),duc_cf(0:imax),ducf(imax)
      real uflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real vflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real wflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real cflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real dfynew_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dfznew_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dfynew_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real dfznew_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc_c (0:imax,jmax/p_row, kmax/p_col)
      real ufluc_s (imax,jmax/p_row, kmax/p_col)
      real cm(imax),cr(imax),qm(imax),qr(imax),Gm(imax),Grr(imax),uc(imax),cm_c(imax),dcc(imax),dc(0:imax)
      real rm(imax),rr(imax),cmf(imax),crf(imax),ucf(imax),cm_cf(imax)
      real em(imax),er(imax),Pam(imax),Pem(imax),Par(imax),Per(imax),Pec(imax),Pac(imax),kGr(imax)
      real kc(imax),ec(imax),wc(imax),Gw(imax),fxc(imax),qint(imax)
      real wmvd(imax),ys(imax),pmtmp(imax),dum(imax),dum5(imax,8),cint(imax),dccc(0:imax),ddc(imax)
      real fxm(imax),fym(imax),fzm(imax),fxr(imax),fyr(imax),fzr(imax)
      real fxmtmp(imax),fymtmp(imax),fzmtmp(imax)
      real qmtmp(imax),Gmtmp(imax),emtmp(imax)
      real dfyc(imax),dfym(imax),dfxm(imax),dfxm_c(0:imax)
      real dfzc(imax),dfzm(imax),subtr(imax),dfxc_c(0:imax),dfxc(imax)
      real interp_c(0:imax,jmax/p_row,kmax/p_col)
      real interp_s(imax,jmax/p_row,kmax/p_col)
      real deriv_s(imax,jmax/p_row,kmax/p_col)
      real epsy(imax),epsz(imax),epsx(imax),ducc(imax),dvcc(imax),dwcc(imax)
      real Rwcz(0:imx,kmax),Rwcy(0:imx,jmax)
      real RGcz(0:imx,kmax),RGcy(0:imx,jmax)
      real Rwcz1(0:imx,kmax),Rwcy1(0:imx,jmax)
      real RGcz1(0:imx,kmax),RGcy1(0:imx,jmax)
      real Rwwz(0:imx,kmax),Rwwy(0:imx,jmax)
      real RGGz(0:imx,kmax),RGGy(0:imx,jmax)
      real Rwwz1(0:imx,kmax),Rwwy1(0:imx,jmax)
      real RGGz1(0:imx,kmax),RGGy1(0:imx,jmax)
      real Rccz(0:imx,kmax),Rccy(0:imx,jmax)
      real Rccz1(0:imx,kmax),Rccy1(0:imx,jmax)
      
      real epsm(imax,8)
      real prodm(imax,8)
      real diffm(imax,8)
      real pdm(imax,8)
      real tdm(imax,8)
      real psm(imax,8)
      real buom(imax,8)
      real radm(imax,8)
      real rdm(imax,8)
      real rem(imax,8)
      real specbudy(0:imx,0:jmax/2,14)
      real specnormy(0:imx,14)
      real specbudz(0:imx,0:kmax/2,14)
      real specnormz(0:imx,14)

      real hi_s(imax,4),hi_c(0:imax,4)
      real ft(jmax,2),dft(jmax,2)
      real fz(kmax,2),dfz(kmax,2)

      real Tref1,Tref
      real R11(imax,kmax),R22(imax,jmax)
      real R11_mpi(imax,kmax),R22_mpi(imax,jmax)

      real speczuu(0:imx,0:kmax/2),speczvv(0:imx,0:kmax/2),speczww(0:imx,0:kmax/2)
      real speczuw(0:imx,0:kmax/2),speczpp(0:imx,0:kmax/2),speczqq(0:imx,0:kmax/2)
      real speczcc(0:imx,0:kmax/2),speczuc(0:imx,0:kmax/2),speczwc(0:imx,0:kmax/2)
      real speczqc(0:imx,0:kmax/2),speczec(0:imx,0:kmax/2),speczkk(0:imx,0:kmax/2)
      real speczee(0:imx,0:kmax/2),speczuc2(0:imx,0:kmax/2),speczcc2(0:imx,0:kmax/2)
      real speczfxc(0:imx,0:kmax/2)
      real speczGG(0:imx,0:kmax/2)
      real specyuu(0:imx,0:jmax/2),specyvv(0:imx,0:jmax/2),specyww(0:imx,0:jmax/2)
      real specyuw(0:imx,0:jmax/2),specypp(0:imx,0:jmax/2),specyqq(0:imx,0:jmax/2)
      real specycc(0:imx,0:jmax/2),specyuc(0:imx,0:jmax/2),specywc(0:imx,0:jmax/2)
      real specyqc(0:imx,0:jmax/2),specyec(0:imx,0:jmax/2),specykk(0:imx,0:jmax/2)
      real specyep1(0:imx,0:jmax/2),specyep2(0:imx,0:jmax/2)
      real specyee(0:imx,0:jmax/2)
      real specyGG(0:imx,0:jmax/2)
      real specyfxc(0:imx,0:kmax/2)
      real crmsy(imax)
      real ecrmsy(imax)
      real qcrmsy(imax)
      real wrmsy(imax)
      real ucrmsy(imax)
      real uwrmsy(imax)
      real wcrmsy(imax)
      real fxcrmsy(imax)
      real energym(imax,4)
      real tempvar(imax,8)
      real turbflu(imax,15)

      real ucQ(imax,4)
      character*5 cha,cha2
      character*1 yes,ycor,yvort
      integer  ire,nfiles,dfiles
      call mpi_init(ierr)
      call decomp_2d_init(imax+2,jmax,kmax,p_row,p_col)

      stime2=MPI_WTIME()

      Rin=0.
      Rout=2.0
      Lz = 16*atan(1.0)
      speczuu=0
      speczvv=0
      speczww=0
      speczuw=0
      speczpp=0
      speczec=0
      speczcc=0
      speczqq=0
      speczqc=0
      speczuc=0
      speczwc=0

      um = 0
      vm = 0
      wm = 0
      rm = 0
      umf = 0
      vmf = 0
      wmf = 0
      pm = 0
      cm = 0
      cmf = 0
      em= 0
      qm = 0
      vcm = 0
      km = 0
      Gm = 0
      fxm= 0
      fym= 0
      fzm= 0
      cr = 0
      ur = 0
      vr = 0
      wr = 0
      rr = 0
      crf = 0
      urf = 0
      vrf = 0
      wrf = 0
      qr = 0
      krr = 0
      Grr = 0
      prr = 0
      ccc = 0
      er = 0
      Pam = 0
      Pem = 0
      Par = 0
      Per = 0
      Pac = 0
      Pec = 0
      subtr = 0      

      tempvar = 0
      turbflu = 0
      energym = 0
      prodm = 0
      epsm = 0
      diffm = 0
      pdm = 0
      tdm =0
      pwm=0
      psm=0
      boum=0
      radm=0
      rdm=0
      rem=0
      specbudy=0
      specnormy=0
      specbudz=0
      specnormz=0
 
      dfyc=0
      dfzc=0

      ucQ = 0
      w3=0
      w4=0
      u3=0
      u4=0
      ctmp=0
      ys= 0
      R11 = 0
      R22 = 0
      rho11 = 0
      rho22 = 0 

      Twt=0.
      Twb=1.
      Gr=0.

      open(5,file='IN')
      read(5,*)       nfiles,iskip,Re,Pr,Pl,Tplus
            close(5)
      if (nrank.eq.0) then
      write(*,222)  nfiles,iskip,Re,Pr,Pl,Tplus
222   format(2I10,' Reynolds= ',F16.3,' Prandtl= ',F16.3,' Planck= ',F16.3,' Lz= ',F16.3, ' Tc= ',F16.3)
      endif

      re = int (RE)
      call cnvstr(ire,cha)
      call mkgrid(Lz,1,nrank)
!      call mk_bin
      do ifiles=1,nfiles
      if (nrank.eq.0) write(6,*) 'reading file number ', ifiles
       call loadd(0,ii,ifiles+iskip)
       call yzderiv(fynew,dfynew_dt,dfynew_dz)
       call yzderiv(fznew,dfznew_dt,dfznew_dz)
        www=0
        ccc=0
        kkk=0
        do k=1,kmax/p_col
         do j=1,jmax/p_row
           do i=0,imax
            ! REYNOLDS AVERAGE
            um(i)=um(i)+unew(i,j,k)/(jmax*kmax)
            ! FAVRE AVERAGE
            umf(i)=umf(i) + phirnew(i,j,k)/(jmax*kmax)
           enddo
           do i=1,imax
             ! REYNOLDS AVERAGES
             vm(i)=vm(i)+vnew(i,j,k)/(jmax*kmax)
             wm(i)=wm(i)+wnew(i,j,k)/(jmax*kmax)
             rm(i)=rm(i)+rhonew(i,j,k)/(jmax*kmax)
             pm(i)=pm(i)+   p(i,j,k)/(jmax*kmax)
             cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax)
             em(i)=em(i)+enew(i,j,k)/(jmax*kmax)
             qm(i)=qm(i)+qrad(i,j,k)/(jmax*kmax)
             vcm(i)=vcm(i)+vcor(i,j,k)/(jmax*kmax)
             fxm(i)=fxm(i)+fxnew(i,j,k)/(jmax*kmax)
             fym(i)=fym(i)+fynew(i,j,k)/(jmax*kmax)
             fzm(i)=fzm(i)+fznew(i,j,k)/(jmax*kmax)
             dfym(i)=dfym(i)+dfynew_dt(i,j,k)/(jmax*kmax)
             dfzm(i)=dfzm(i)+dfznew_dz(i,j,k)/(jmax*kmax)
             Gm(i)=Gm(i)+Gnew(i,j,k)/(jmax*kmax)
             km(i)=km(i)+knew(i,j,k)/(jmax*kmax)
             Pem(i)=Pem(i)+knew(i,j,k)*enew(i,j,k)/(jmax*kmax)
             Pam(i)=Pam(i)+knew(i,j,k)*Gnew(i,j,k)/(jmax*kmax)
             ccc(i)=ccc(i)+cnew(i,j,k)*wnew(i,j,k)/(jmax*kmax)
             www(i)=www(i)+wnew(i,j,k)/(jmax*kmax)
             kkk(i)=kkk(i)+knew(i,j,k)/(jmax*kmax)
             ! FAVRE AVERAGES
             vmf(i)=vmf(i)+phitnew(i,j,k)/(jmax*kmax)
             wmf(i)=wmf(i)+phiznew(i,j,k)/(jmax*kmax)
             cmf(i)=cmf(i)+  renew(i,j,k)/(jmax*kmax)             
            enddo
          enddo
        enddo

      call mpi_allreduce(www,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      www = mpi_sm

      call mpi_allreduce(ccc,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ccc = mpi_sm

      call mpi_allreduce(kkk,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kkk = mpi_sm
 
      bulk1 = 0
      cbulk = 0
      opthk = 0 
       do i=1,imax
        bulk1 = bulk1 + (Ru(i)-Ru(i-1))/2.*www(i)
        cbulk = cbulk + (Ru(i)-Ru(i-1))/2.*ccc(i)
        opthk = opthk + (Ru(i)-Ru(i-1))/2.*kkk(i)
       enddo
       cbulk=cbulk/bulk1
       if (nrank.eq.0) write(6,21) ifiles,bulk1,cbulk,opthk 
21     format( 'file= ',I5,' bulk= ',F10.6,' cbulk= ',F10.6,' opthk= ',F10.6)
      enddo

!------------------------------------ END READING FILES


      call mpi_allreduce(um,mpi_cm,imax+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      um = mpi_cm /(nfiles)
      call mpi_allreduce(vm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vm = mpi_sm /(nfiles)
      call mpi_allreduce(wm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wm = mpi_sm /(nfiles)
      call mpi_allreduce(cm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      cm = mpi_sm /(nfiles)
      call mpi_allreduce(rm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      rm = mpi_sm /(nfiles)
      call mpi_allreduce(pm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      pm = mpi_sm /(nfiles)
      call mpi_allreduce(qm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qm = mpi_sm/(nfiles)    
      call mpi_allreduce(vcm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vcm = mpi_sm/(nfiles)    
      call mpi_allreduce(fxm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      fxm = mpi_sm/(nfiles)    
      call mpi_allreduce(fym,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      fym = mpi_sm/(nfiles)    
      call mpi_allreduce(fzm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      fzm = mpi_sm/(nfiles)    
      call mpi_allreduce(dfym,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      dfym = mpi_sm/(nfiles)    
      call mpi_allreduce(dfzm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      dfzm = mpi_sm/(nfiles)    
      call mpi_allreduce(km,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      km = mpi_sm/(nfiles)    
      call mpi_allreduce(Gm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gm = mpi_sm/(nfiles)    
      call mpi_allreduce(em,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      em = mpi_sm/(nfiles)
      call mpi_allreduce(Pam,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Pam = mpi_sm/(nfiles)
      call mpi_allreduce(Pem,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Pem = mpi_sm/(nfiles)


      call mpi_allreduce(umf,mpi_cm,imax+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      umf = mpi_cm /(nfiles)
      call mpi_allreduce(vmf,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vmf = mpi_sm /(nfiles * rm)
      call mpi_allreduce(wmf,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wmf = mpi_sm /(nfiles * rm)
      call mpi_allreduce(cmf,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      cmf = mpi_sm /(nfiles * rm)
     
      mpi_cm = 0
      call inter_s_c_m(imax,1,rm,mpi_cm,dr)
      umf = umf/mpi_cm

      call der1w_s_c6_m(imax,1,wm,dw,dr)
      call der1c_s_c6_m(imax,1,cm,dc,dr,Twb,Twt)
      dw=dw*mr_c
      dc=dc*mr_c
      call inter_c_s_m(imax,1,dw,dww,dr)
      call inter_c_s_m(imax,1,dc,dcc,dr)
 
      call deriv_s_c_m(imax,1,dcc,dccc,dr)
      dccc=dccc*mr_c
      call inter_c_s_m(imax,1,dccc,ddc,dr)


      do ifiles=1,nfiles
      if (nrank.eq.0) write(6,*) 'reading file number ', ifiles
      call loadd(0,ii,ifiles+iskip)
      pmtmp=0
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         pmtmp(i)=pmtmp(i)+   p(i,j,k)
        enddo
       enddo
      enddo
       pmtmp=pmtmp/(jmax*kmax)

       call mpi_allreduce(pmtmp,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
       pmtmp = mpi_sm

       do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=0,imax
         ! REYNOLDS FLUCTUATIONS
         ufluc(i,j,k) = unew(i,j,k) - um(i)
         uMean(i,j,k) = um(i)
         ! FAVRE FLUCTUATIONS
         uflucf(i,j,k)= unew(i,j,k) - umf(i)
        enddo
        do i=1,imax
         vMean(i,j,k) = vm(i)
         wMean(i,j,k) = wm(i)
         ! REYNOLDS FLUCTUATIONS
         vfluc(i,j,k) = vnew(i,j,k) - vm(i)
         wfluc(i,j,k) = wnew(i,j,k) - wm(i)
         rfluc(i,j,k) = rhonew(i,j,k) - rm(i)
         pfluc(i,j,k) = p(i,j,k) - pmtmp(i)
         cfluc(i,j,k) = cnew(i,j,k) - cm(i)
         qfluc(i,j,k) = qrad(i,j,k) - qm(i)
         Gfluc(i,j,k) = Gnew(i,j,k) - Gm(i)
         Pafluc(i,j,k) = knew(i,j,k)*Gnew(i,j,k) - Pam(i)
         Pefluc(i,j,k) = knew(i,j,k)*enew(i,j,k) - Pem(i)
         efluc(i,j,k) = enew(i,j,k) - em(i)
         kfluc(i,j,k) = knew(i,j,k) - km(i)
         fxfluc(i,j,k) = fxnew(i,j,k) - fxm(i)
         fyfluc(i,j,k) = fynew(i,j,k) - fym(i)
         fzfluc(i,j,k) = fznew(i,j,k) - fzm(i)
         ! FAVRE FLUCTUATIONS
         vflucf(i,j,k) = vnew(i,j,k) - vmf(i)
         wflucf(i,j,k) = wnew(i,j,k) - wmf(i)
         cflucf(i,j,k) = cnew(i,j,k) - cmf(i)
        enddo
       enddo
      enddo
      do i=0,imax
      ufluc_c(i,:,:)=ufluc(i,:,:)
      enddo

      ufluc_tmp=0
      dcfluc=0
      dc_dt=0
      dc_dz=0

      call inter_c_s_m(imax,jmax*kmax/px,ufluc_c,ufluc_s,dr)
      do i=1,imax
       ufluc_tmp(i,:,:)=ufluc_s(i,:,:)  
      enddo 

      call yzderiv(cfluc,dc_dt,dc_dz)
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dcfluc(i,:,:)=deriv_s(i,:,:)
      enddo
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)*ufluc_s(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       ducfluc(i,:,:)=deriv_s(i,:,:)
      enddo
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)*vfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dvcfluc(i,:,:)=deriv_s(i,:,:)
      enddo
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)*wfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dwcfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      call meanenergy(nfiles,rm,umf,cm,cmf,uflucf,cflucf,fxm,energym)
      call tempvarbal(nfiles,rm,umf,cm,cmf,uflucf,cfluc,cflucf,qm,qfluc,tempvar)
      call turbflubal(nfiles,rm,um,umf,cm,cmf,qm,pmtmp,ufluc,uflucf,cfluc,cflucf,qfluc,
     +                pfluc,vflucf,wflucf,vfluc,wfluc,turbflu)

      call bal(um,vm,wm,pmtmp,cm,ufluc,vfluc,wfluc,pfluc,cfluc,ufluc_s,qfluc,fxfluc,fyfluc,fzfluc,
     +         epsm,prodm,diffm,tdm,pdm,psm,buom,radm,rdm,rem,dfyc,dfzc,nfiles,nrank)

      call hos(um,vm,wm,pmtmp,cm,em,qm,Gm,km,Pam,Pem,Par,Per,Pec,Pac,kGr,ufluc_s,ur,vr,wr,cr,er,qr,prr,Grr,krr,
     +        uu,vv,ww,uw,uc,wc,kc,ec,fxr,fyr,fzr,fxm,fym,fzm,fxc,dcfluc,dc_dt,dc_dz,epsx,epsy,epsz,Gw,
     +        ducfluc,dvcfluc,dwcfluc,ducc,dvcc,dwcc,umf,vmf,wmf,cmf,urf,vrf,wrf,crf,nfiles)

      if(nrank.eq.0) write(*,*) 'entering specbalz' 
      stime=MPI_WTIME()
      call spectrazcalc(dcc,ufluc,ufluc_tmp,vfluc,wfluc,cfluc,pfluc,qfluc,efluc,Pefluc,Pafluc,kfluc,fxfluc,fyfluc,fzfluc,
     +         Gfluc,speczuu,speczGG,
     +         speczvv,speczww,speczuw,speczuc,speczwc,speczqq,speczqc,speczcc,speczec,speczee,speczkk,speczfxc,specbudz,
     +               nfiles)

      if(nrank.eq.0) write(*,*) 'entering specbaly' 

      call spectraycalc(dcc,ufluc,ufluc_tmp,vfluc,wfluc,cfluc,pfluc,qfluc,efluc,Pefluc,Pafluc,kfluc,fxfluc,fyfluc,fzfluc,
     +         Gfluc,specyuu,specyGG,
     +         specyvv,specyww,specyuw,specyuc,specywc,specyqq,specyqc,specycc,specyec,specyee,specykk,specyfxc,specbudy,
     +               nfiles)

      call corry(wfluc,cfluc,Rwcy1,nfiles)
      call corry(Gfluc,cfluc,RGcy1,nfiles)
      call corry(cfluc,cfluc,Rccy1,nfiles)
      call corry(Gfluc,Gfluc,RGGy1,nfiles)
      call corry(wfluc,wfluc,Rwwy1,nfiles)
      call corrz(wfluc,cfluc,Rwcz1,nfiles)
      call corrz(Gfluc,cfluc,RGcz1,nfiles)
      call corrz(cfluc,cfluc,Rccz1,nfiles)
      call corrz(Gfluc,Gfluc,RGGz1,nfiles)
      call corrz(wfluc,wfluc,Rwwz1,nfiles)

      call quad_analysis(ucQ,ufluc_tmp,cfluc,nfiles)

      do i=0,imx
       do j=1,jmax
        Rwcy(i,j)=Rwcy(i,j)+Rwcy1(i,j)/nfiles
        RGcy(i,j)=RGcy(i,j)+RGcy1(i,j)/nfiles
        Rccy(i,j)=Rccy(i,j)+Rccy1(i,j)/nfiles
        RGGy(i,j)=RGGy(i,j)+RGGy1(i,j)/nfiles
        Rwwy(i,j)=Rwwy(i,j)+Rwwy1(i,j)/nfiles
       enddo
       do k=1,kmax
        Rwcz(i,k)=Rwcz(i,k)+Rwcz1(i,k)/nfiles
        RGcz(i,k)=RGcz(i,k)+RGcz1(i,k)/nfiles
        Rccz(i,k)=Rccz(i,k)+Rccz1(i,k)/nfiles
        RGGz(i,k)=RGGz(i,k)+RGGz1(i,k)/nfiles
        Rwwz(i,k)=Rwwz(i,k)+Rwwz1(i,k)/nfiles
       enddo
      enddo
 

      if(nrank.eq.0) write(*,*) 'time in specbal: ',MPI_WTIME()-stime 
     
    
      do i=1,imax
       subtr(i)=subtr(i)+(uc(i)-dcc(i)/(Pr*Re)+fxmtmp(i)/(Re*Pr*Pl))/nfiles
      enddo

      enddo


      call parsum(ur,vr,wr,cr,prr,u3,ww)
      call parsum(urf,vrf,wrf,crf,dum,dum,dum)
      call parsum(v3,w3,u4,v4,w4,uu,vv)
      call parsum(uw,uc,wc,qr,Grr,krr,kGr)
      call parsum(er,Par,Per,Pec,Pac,kc,ec)
      call parsum(fxr,fyr,fzr,dfyc,dfzc,fxc,dum)
      call parsum(epsx,epsy,epsz,ducc,dum,dum,dum)

      call parsum5(epsm,prodm,diffm,tdm,pdm,psm,buom)
      call parsum5(radm,rdm,tempvar,dum5,dum5,dum5,dum5)

      call parsum4(energym,ucQ)
      call parsum15(turbflu)

      call deriv_s_c_m(imax,1,uc,duc_c,dr)
      duc_c=duc_c*mr_c
      call inter_c_s_m(imax,1,duc_c,duc,dr) 

      call deriv_s_c_m(imax,1,fxm,dfxm_c,dr)
      dfxm_c=dfxm_c*mr_c
      call inter_c_s_m(imax,1,dfxm_c,dfxm,dr) 

      call deriv_s_c_m(imax,1,fxc,dfxc_c,dr)
      dfxc_c=dfxc_c*mr_c
      call inter_c_s_m(imax,1,dfxc_c,dfxc,dr)

      if(nrank.eq.0) write(*,*) 'calculating normalization factors'
      
      do rs=1,14
       do i=1,imax
        do j=0,jmax/2
         specnormy(i,rs)=specnormy(i,rs)+specbudy(i,j,rs)
        enddo
        do k=0,kmax/2
         specnormz(i,rs)=specnormz(i,rs)+specbudz(i,k,rs)
        enddo
       enddo
      enddo
     
      qint(1)=-qm(1)*(ru(i)-ru(i-1))
      do i=2,imax
       qint(i)=qint(i-1)-qm(i)*(ru(i)-ru(i-1))
      enddo

      if (nrank.eq.0)  write(6,*) 'done with : reading the data'
      
      if (nrank.eq.0) then

      open(31,file='Results/specbudy')
      do i=1,imax
       do j=0,jmax/2
       write(31,'(25E16.5)') rp(i),j*2./1.5,specbudy(i,j,1),specbudy(i,j,2),specbudy(i,j,3),specbudy(i,j,4),
     +                    specbudy(i,j,5),specbudy(i,j,6),specbudy(i,j,7),specbudy(i,j,8),specbudy(i,j,9),
     +                    specbudy(i,j,10),specbudy(i,j,11),specbudy(i,j,12),specbudy(i,j,13),specbudy(i,j,14)
       enddo
      enddo
      close(31)

      open(31,file='Results/specbudz')
      do i=1,imax
       do k=0,kmax/2
       write(31,'(25E16.5)') rp(i),k*2./4.,specbudz(i,k,1),specbudz(i,k,2),specbudz(i,k,3),specbudz(i,k,4),
     +                    specbudz(i,k,5),specbudz(i,k,6),specbudz(i,k,7),specbudz(i,k,8),specbudz(i,k,9),
     +                    specbudz(i,k,10),specbudz(i,k,11),specbudz(i,k,12),specbudz(i,k,13),specbudz(i,k,14)
       enddo
      enddo
      close(31)

      open(31,file='Results/specnormy')
      do i=1,imax
       write(31,'(25E16.5)') rp(i),specnormy(i,1),specnormy(i,2),specnormy(i,3),specnormy(i,4),specnormy(i,5),specnormy(i,6),
     +                    specnormy(i,7),specnormy(i,8),specnormy(i,9),specnormy(i,10),specnormy(i,11),specnormy(i,12)
     +                   ,specnormy(i,13),specnormy(i,14) 
      enddo
      close(31)

      open(31,file='Results/specnormz')
      do i=1,imax
       write(31,'(25E16.5)') rp(i),specnormz(i,1),specnormz(i,2),specnormz(i,3),specnormz(i,4),specnormz(i,5),specnormz(i,6),
     +                     specnormz(i,7),specnormz(i,8),specnormz(i,9),specnormz(i,10),specnormz(i,11),specnormz(i,12)
     +                    ,specnormz(i,13),specnormz(i,14) 
      enddo
      close(31)

      open(31,file='Results/kbud')
      do i=1,imax
       write(31,'(14E16.5)') rp(i),(ur(i)+vr(i)+wr(i))/2.,epsm(i,1),prodm(i,1),psm(i,1),pdm(i,1),tdm(i,1),diffm(i,1),radm(i,1)
      enddo
      close(31)

      open(32,file='Results/ubud')
      do i=1,imax
       write(33,'(14E16.5)') rp(i),epsm(i,3),prodm(i,3),psm(i,3),pdm(i,3),tdm(i,3),diffm(i,3),radm(i,3)
      enddo
      close(33)


      open(34,file='Results/wbud')
      do i=1,imax
       write(34,'(14E16.5)') rp(i),epsm(i,4),prodm(i,4),psm(i,4),pdm(i,4),tdm(i,4),diffm(i,4),radm(i,4)
      enddo
      close(34)

      open(35,file='Results/uwbud')
      do i=1,imax
       write(35,'(14E16.5)') rp(i),epsm(i,5),prodm(i,5),psm(i,5),pdm(i,5),tdm(i,5),diffm(i,5),radm(i,5)
      enddo
      close(35)

      open(35,file='Results/cbud')
      do i=1,imax
       write(35,'(14E16.5)') rp(i),cr(i),epsm(i,6),prodm(i,6),psm(i,6),pdm(i,6),tdm(i,6),diffm(i,6),radm(i,6),
     +                       rdm(i,6),rem(i,6),Pac(i)/(Re*Pr*Pl)*2,Pec(i)*2/(Re*Pr*Pl)
      enddo
      close(35)

      open(35,file='Results/tempvar')
      do i=1,imax
       write(35,'(14E16.5)') rp(i),crf(i),tempvar(i,1),tempvar(i,2),
     +                       tempvar(i,3),tempvar(i,4),tempvar(i,5),
     +                       tempvar(i,6),tempvar(i,7),tempvar(i,8)
      enddo
      close(35)
    
      open(35,file='Results/turbflu')
      do i=1,imax
       write(35,'(20E16.5)') rp(i),uc(i), turbflu(i,1),turbflu(i,2),
     +                       turbflu(i,3),turbflu(i,4),turbflu(i,5),
     +                       turbflu(i,6),turbflu(i,7),turbflu(i,8),
     +                       turbflu(i,9),turbflu(i,10),turbflu(i,11),
     +                       turbflu(i,12),turbflu(i,13),turbflu(i,14),
     +                       turbflu(i,15)
      enddo
      close(35)


      open(35,file='Results/wcbud')
      do i=1,imax
       write(35,'(14E16.5)') rp(i),epsm(i,7),prodm(i,7),psm(i,7),pdm(i,7),tdm(i,7),diffm(i,7),radm(i,7),
     +                       rdm(i,7),rem(i,7)
      enddo
      close(35)

      open(35,file='Results/ucbud')
      do i=1,imax
       write(35,'(14E16.5)') rp(i),epsm(i,8),prodm(i,8),psm(i,8),pdm(i,8),tdm(i,8),diffm(i,8),radm(i,8),
     +                       rdm(i,8),rem(i,8)
      enddo
      close(35)

      open(36,file='Results/meanenergy')
      do i=1,imax
       write(36,'(24E16.5)') rp(i),energym(i,1),energym(i,2),energym(i,3),energym(i,4) 
      end do
      close(36)

      open(36,file='Results/resid')
      do i=1,imax
       write(36,'(24E16.5)') rp(i),sqrt(ur(i)),sqrt(vr(i)),sqrt(wr(i)),sqrt(cr(i)),sqrt(prr(i)),sqrt(qr(i)),sqrt(er(i))
     +                       ,sqrt(Grr(i)),sqrt(krr(i)),sqrt(Par(i)),sqrt(Per(i)),sqrt(fxr(i)),sqrt(fyr(i)),sqrt(fzr(i)),
     +                        kGr(i),kc(i),ec(i)
      end do
      close(36)

      open(36,file='Results/residf')
      do i=1,imax
       write(36,'(24E16.5)') rp(i),sqrt(urf(i)),sqrt(vrf(i)),sqrt(wrf(i)),sqrt(crf(i))
      end do
      close(36)

      open(36,file='Results/radiation')
      do i=1,imax
       write(36,'(20E16.5)') rp(i),sqrt(cr(i)),sqrt(qr(i)),sqrt(er(i)),sqrt(Grr(i)),sqrt(krr(i)),
     +                       sqrt(Par(i)),sqrt(Per(i)),kGr(i),Pec(i),Pac(i),kc(i),ec(i),em(i),cm(i),km(i),Gm(i)
      end do
      close(36)
  

      open(36,file='Results/mean')
      do i=1,imax
       write(36,'(21E16.5)') rp(i),um(i),vm(i),wm(i),cm(i),pm(i),qm(i),vcm(i),em(i),km(i),Gm(i),Pam(i),Pem(i),rm(i),
     +                       uu(i),vv(i),ww(i),uw(i)
      end do
      close(36)

      open(36,file='Results/meanf')
      do i=1,imax
       write(36,'(21E16.5)') rp(i),umf(i),vmf(i),wmf(i),cmf(i)
      end do
      close(36)

      open(36,file='Results/stress')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),uw(i),dww(i)/Re
      end do
      close(36)


      open(36,file='Results/quadrant')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),ucQ(i,1),ucQ(i,2),ucQ(i,3),ucQ(i,4)
      end do
      close(36)


      open(36,file='Results/flux')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),uc(i),dcc(i)/(Pr*Re),fxm(i)/(Re*Pr*Pl),qint(i)/(Re*Pr*Pl),wc(i),Gw(i)
      end do
      close(36)

      open(36,file='Results/fluxrad')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),fxm(i),fym(i),fzm(i),dfym(i),dfzm(i),dfxm(i),dfyc(i),dfzc(i),qm(i)
      end do
      close(36)

      open(36,file='Results/model')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),uc(i),uw(i),dcc(i),dww(i),fxm(i)
      end do
      close(36)

      open(36,file='Results/corrz')
      do i=0,imx
       do k=1,kmax
        write(36,'(25E16.5)') rp(i),k*1.,Rwcz(i,k),RGcz(i,k),Rccz(i,k),RGGz(i,k),Rwwz(i,k)
       enddo
      end do
      close(36)

      open(36,file='Results/corry')
      do i=0,imx
       do j=1,jmax
        write(36,'(25E16.5)') rp(i),j*1.,Rwcy(i,j),RGcy(i,j),Rccy(i,j),RGGy(i,j),Rwwy(i,j)
       enddo
      end do
      close(36)

      open(36,file='Results/fluxbal')
      do i=1,imax
       write(36,'(14E16.5)') rp(i),duc(i),ddc(i)/(Pr*Re),qm(i)/(Re*Pr*Pl)
      end do
      close(36)

      open(36,file='Results/spectrazv')
      do k=1,kmax/2-1
        write(36,'(10E20.6)') (k-1)*2./4.,speczuu(imx-11,k),speczvv(imx-11,k),speczww(imx-11,k),speczuw(imx-11,k)
      end do
      close(36)

      open(36,file='Results/spectrazc')
      do k=1,kmax/2-1
        write(36,'(10E20.6)') (k-1)*2./4.,speczuc(imx-11,k),speczwc(imx-11,k),speczcc(imx-11,k),speczqc(imx-11,k)
      end do
      close(36)

      open(36,file='Results/spectrazo')
      do k=1,kmax/2-1
        write(36,'(10E20.6)') (k-1)*2./4.,speczqq(imx-11,k),speczpp(imx-11,k),speczec(imx-11,k),speczkk(imx-11,k),
     +                                    speczee(imx-11,k)
      end do
      close(36)

      open(36,file='Results/spectrayv')
      do j=1,jmax/2-1
        write(36,'(10E20.6)') (j-1)*2./1.5,specyuu(imx-11,j),specyvv(imx-11,j),specyww(imx-11,j),specyuw(imx-11,j)
      end do
      close(36)

      open(36,file='Results/spectrayc')
      do j=1,jmax/2-1
        write(36,'(10E20.6)') (j-1)*2./1.5,specyuc(imx-11,j),specywc(imx-11,j),specycc(imx-11,j),specyqc(imx-11,j)
      end do
      close(36)

      open(36,file='Results/spectrayo')
      do j=1,jmax/2-1
        write(36,'(10E20.6)') (j-1)*2./1.5,specyqq(imx-11,j),specypp(imx-11,j),specyec(imx-11,j),specykk(imx-11,j),
     +                                     specyee(imx-11,j)
      end do
      close(36)

      open(36,file='Results/spectraycomp')
      do i=1,imax
       do j=1,jmax/2-1
        write(36,'(10E20.6)') rp(i),(j)*2./1.5,specyqq(i,j),specypp(i,j),specyec(i,j),specykk(i,j),
     +                                     specyee(i,j),specycc(i,j),specyqc(i,j),specyGG(i,j)
       enddo
      enddo
      close(36)

      open(36,file='Results/spectrazcomp')
      do i=1,imax
       do k=1,kmax/2-1
        write(36,'(10E20.6)') rp(i),(k)*2./4.,speczqq(i,k),speczpp(i,k),speczec(i,k),speczkk(i,k),
     +                                     speczee(i,k),speczcc(i,k),speczqc(i,k),speczGG(i,k)
       enddo
      enddo
      close(36)



      write(*,*) 'Time required: ',MPI_WTIME()-stime2

      end if
!-------------------------------------- END PROGRAM

      call decomp_2d_finalize
      call mpi_finalize(ierr)

      end

      subroutine hos(um,vm,wm,pmtmp,cm,em,qm,Gm,km,Pam,Pem,Par,Per,Pec,Pac,kGr,ufluc_s,ur,vr,wr,cr,er,qr,prr,Grr,krr,
     +        uu,vv,ww,uw,uc,wc,kc,ec,fxr,fyr,fzr,fxm,fym,fzm,fxc,dcfluc,dc_dt,dc_dz,epsx,epsy,epsz,Gw,
     +        ducfluc,dvcfluc,dwcfluc,ducc,dvcc,dwcc,umf,vmf,wmf,cmf,urf,vrf,wrf,crf,nfiles)

      implicit none
      include 'par_post.txt'
      include 'common.txt'
      real um(0:imax),vm(imax),wm(imax),pmtmp(imax),prr(imax),qm(imax),qr(imax)
      real umf(0:imax),vmf(imax),wmf(imax),cmf(imax)
      real urf(0:imax),vrf(imax),wrf(imax),crf(imax)
      real em(imax),km(imax),Gm(imax),Grr(imax),krr(imax),er(imax),wc(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),dw(0:imax),dww(imax),w3(imax),w4(imax)
      real u3(imax),u4(imax),v3(imax),v4(imax),bulk1,www(imax),bulkrms
      real Pam(imax),Pem(imax),Par(imax),Per(imax),kGr(imax),Pec(imax),Pac(imax),Gw(imax)
      real cm(imax),cr(imax),cm_c(imax),dcc(imax),dc(0:imax),uc(imax),kc(imax),ec(imax)
      real uu(imax),vv(imax),ww(imax),ufluc_s(imax,jmax/p_row,kmax/p_col),epsx(imax),epsy(imax),epsz(imax)
      real fxr(imax),fyr(imax),fzr(imax),fxm(imax),fym(imax),fzm(imax),fxc(imax),ducc(imax),dvcc(imax),dwcc(imax)
      real dcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real ducfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dvcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dwcfluc(0:imax+1,jmax/p_row, kmax/p_col)

      integer nfiles

         do k=1,kmax/p_col
          do j=1,jmax/p_row
            do i=1,imax
             ! REYNOLDS RMS VALUES
             ur(i)=ur(i)+(unew(i,j,k)-um(i))**2/((jmax)*(kmax)*nfiles)
             vr(i)=vr(i)+(vnew(i,j,k)-vm(i))**2/((jmax)*(kmax)*nfiles)
             wr(i)=wr(i)+(wnew(i,j,k)-wm(i))**2/((jmax)*(kmax)*nfiles)
             cr(i)=cr(i)+(cnew(i,j,k)-cm(i))**2/((jmax)*(kmax)*nfiles)
             er(i)=er(i)+(enew(i,j,k)-em(i))**2/((jmax)*(kmax)*nfiles)
             qr(i)=qr(i)+(qrad(i,j,k)-qm(i))**2/((jmax)*(kmax)*nfiles)
             fxr(i)=fxr(i)+(fxnew(i,j,k)-fxm(i))**2/((jmax)*(kmax)*nfiles)
             fyr(i)=fyr(i)+(fynew(i,j,k)-fym(i))**2/((jmax)*(kmax)*nfiles)
             fzr(i)=fzr(i)+(fznew(i,j,k)-fzm(i))**2/((jmax)*(kmax)*nfiles)
             !uc(i)=uc(i)+(cnew(i,j,k)-cm(i))*ufluc_s(i,j,k)/((jmax)*(kmax)*nfiles)
             epsx(i)=epsx(i)+(fxnew(i,j,k)-fxm(i))*dcfluc(i,j,k)/((jmax)*(kmax)*nfiles)
             epsy(i)=epsy(i)+(fynew(i,j,k)-fym(i))*dc_dt(i,j,k)/((jmax)*(kmax)*nfiles)
             epsz(i)=epsz(i)+(fznew(i,j,k)-fzm(i))*dc_dz(i,j,k)/((jmax)*(kmax)*nfiles)
             wc(i)=wc(i)+(cnew(i,j,k)-cm(i))*(wnew(i,j,k)-wm(i))/((jmax)*(kmax)*nfiles)
             Gw(i)=Gw(i)+(Gnew(i,j,k)-Gm(i))*(wnew(i,j,k)-wm(i))/((jmax)*(kmax)*nfiles)
             fxc(i)=fxc(i)+(cnew(i,j,k)-cm(i))*(fxnew(i,j,k)-fxm(i))/((jmax)*(kmax)*nfiles)
             ducc(i)=ducc(i)+ducfluc(i,j,k)*(cnew(i,j,k)-cm(i))/((jmax)*(kmax)*nfiles)
             dvcc(i)=dvcc(i)+dvcfluc(i,j,k)*(cnew(i,j,k)-cm(i))/((jmax)*(kmax)*nfiles)
             dwcc(i)=dwcc(i)+dwcfluc(i,j,k)*(cnew(i,j,k)-cm(i))/((jmax)*(kmax)*nfiles)

             Grr(i)=Grr(i)+(Gnew(i,j,k)-Gm(i))**2/((jmax*kmax)*nfiles)
             krr(i)=krr(i)+(knew(i,j,k)-km(i))**2/((jmax*kmax)*nfiles)
             kGr(i)=kGr(i)+(knew(i,j,k)-km(i))*(Gnew(i,j,k)-Gm(i))/((jmax*kmax)*nfiles)
             Par(i)=Par(i)+(knew(i,j,k)*Gnew(i,j,k)-Pam(i))**2/((jmax*kmax)*nfiles)
             Per(i)=Per(i)+(knew(i,j,k)*enew(i,j,k)-Pem(i))**2/((jmax*kmax)*nfiles)
     
             Pac(i)=Pac(i)+(knew(i,j,k)*Gnew(i,j,k)-Pam(i))*(cnew(i,j,k)-cm(i))/((jmax*kmax)*nfiles)
             Pec(i)=Pec(i)+(knew(i,j,k)*enew(i,j,k)-Pem(i))*(cnew(i,j,k)-cm(i))/((jmax*kmax)*nfiles)

             kc(i)=kc(i)+(knew(i,j,k)-km(i))*(cnew(i,j,k)-cm(i))/(jmax*kmax*nfiles) 
             ec(i)=ec(i)+(enew(i,j,k)-em(i))*(cnew(i,j,k)-cm(i))/(jmax*kmax*nfiles)
             prr(i)=prr(i)+(p(i,j,k)-pmtmp(i))**2/((jmax)*(kmax)*nfiles)
             w3(i)=w3(i)+(wnew(i,j,k)-wm(i))**3/((jmax)*(kmax)*nfiles)
             w4(i)=w4(i)+(wnew(i,j,k)-wm(i))**4/((jmax)*(kmax)*nfiles)
             v3(i)=v3(i)+(vnew(i,j,k)-vm(i))**3/((jmax)*(kmax)*nfiles)
             v4(i)=v4(i)+(vnew(i,j,k)-vm(i))**4/((jmax)*(kmax)*nfiles)
             u3(i)=u3(i)+(unew(i,j,k)-um(i))**3/((jmax)*(kmax)*nfiles)
             u4(i)=u4(i)+(unew(i,j,k)-um(i))**4/((jmax)*(kmax)*nfiles)
             uw(i)=uw(i)+(wnew(i,j,k)-wm(i))*ufluc_s(i,j,k)/((jmax)*(kmax)*nfiles)
             uu(i)=uu(i)+ufluc_s(i,j,k)*ufluc_s(i,j,k)/((jmax)*(kmax)*nfiles)
             vv(i)=vv(i)+(vnew(i,j,k)-vm(i))**2./((jmax)*(kmax)*nfiles)
             ww(i)=ww(i)+(wnew(i,j,k)-wm(i))**2./((jmax)*(kmax)*nfiles)

             ! FAVRE RMS VALUES
             urf(i)=urf(i)+(unew(i,j,k)-umf(i))**2/((jmax)*(kmax)*nfiles)
             vrf(i)=vrf(i)+(vnew(i,j,k)-vmf(i))**2/((jmax)*(kmax)*nfiles)
             wrf(i)=wrf(i)+(wnew(i,j,k)-wmf(i))**2/((jmax)*(kmax)*nfiles)
             crf(i)=crf(i)+(cnew(i,j,k)-cmf(i))**2/((jmax)*(kmax)*nfiles)
             wc(i)=wc(i)+rhonew(i,j,k)*(cnew(i,j,k)-cm(i))*(wnew(i,j,k)-wm(i))/((jmax)*(kmax)*nfiles)
             uc(i)=uc(i)+rhonew(i,j,k)*(cnew(i,j,k)-cm(i))*(unew(i,j,k)+unew(i-1,j,k)-um(i)-um(i-1))*.5/((jmax)*(kmax)*nfiles)

            enddo
          enddo
        enddo

      end

      subroutine bal(um,vm,wm,pm,cm,ufluc,vfluc,wfluc,pfluc,cfluc,ufluc_s,qfluc,fxfluc,fyfluc,fzfluc,
     +         epsm,prodm,diffm,tdm,pdm,psm,buom,radm,rdm,rem,dfyc,dfzc,nf,nrank)
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'

      real um(0:imax),vm(imax),wm(imax),pm(imax),pm_c(0:imax),dpm_s(imax),cm(imax)
      real eps(0:imax+1,jmax/p_row, kmax/p_col,8)
      real epsrad(0:imax+1,jmax/p_row, kmax/p_col,8)
      real prod(0:imax+1,jmax/p_row, kmax/p_col,8)
      real ps(0:imax+1,jmax/p_row,kmax/p_col,8)
      real ufluc(0:imax+1,jmax/p_row, kmax/p_col)
      real vfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real pfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real cfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real qfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fxfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fyfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fzfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc_s (imax    ,jmax/p_row, kmax/p_col)
      real tke_s(imax,jmax/p_row, kmax/p_col,8)
      real pac_s(imax,jmax/p_row,kmax/p_col,8)
      real rd_s(imax,jmax/p_row,kmax/p_col,8)
      real buo(imax,jmax/p_row,kmax/p_col,8)
      real rad1(imax,jmax/p_row,kmax/p_col,8)
      real diff(imax,jmax/p_row,kmax/p_col,8)
      real diff_s(imax,jmax/p_row,kmax/p_col,8)
      real difff(0:imax+1,jmax/p_row,kmax/p_col,8)
      real difff2(0:imax+1,jmax/p_row,kmax/p_col,8)
      real diff_c(0:imax,jmax/p_row,kmax/p_col,8)
      real dtke_s(imax,jmax/p_row, kmax/p_col,8)
      real dpac_s(imax,jmax/p_row,kmax/p_col,8)
      real drd_s(imax,jmax/p_row,kmax/p_col,8)
      real tke_c(0:imax,jmax/p_row, kmax/p_col,8)
      real pac_c(0:imax,jmax/p_row,kmax/p_col,8)
      real rd_c(0:imax,jmax/p_row,kmax/p_col,8)
      real hiums_s(imax, jmax/p_row, kmax/p_col,8),hiums_c(0:imax, jmax/p_row, kmax/p_col,8)
      real test_s(imax),test_c(0:imax)
      real dfyc_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dfzc_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dfyc(imax)
      real dfzc(imax)

      real epsm(imax,8)
      real prodm(imax,8)
      real diffm(imax,8)
      real tdm(imax,8)
      real pdm(imax,8)
      real psm(imax,8)
      real buom(imax,8)
      real radm(imax,8)
      real rdm(imax,8)
      real rem(imax,8)
      integer nf,rs,nrank,ier



      real utmp(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmp(0:imax+1,jmax/p_row,kmax/p_col)

      epsrad = 0.
      rd_s = 0. 
      rad1(:,:,:,:) = 0
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax

          pac_s(i,j,k,3) = 0
          pac_s(i,j,k,4) = 0
          pac_s(i,j,k,1) = pfluc(i,j,k)*ufluc_s(i,j,k)
          pac_s(i,j,k,2) = 2*pfluc(i,j,k)*ufluc_s(i,j,k)
          pac_s(i,j,k,5) = pfluc(i,j,k)*wfluc(i,j,k)
          pac_s(i,j,k,6) = 0
          pac_s(i,j,k,7) = 0
          pac_s(i,j,k,8) = pfluc(i,j,k)*cfluc(i,j,k)

          tke_s(i,j,k,2) = (ufluc_s(i,j,k)**2)*ufluc_s(i,j,k)
          tke_s(i,j,k,3) = (vfluc(i,j,k)**2)*ufluc_s(i,j,k)
          tke_s(i,j,k,4) = (wfluc(i,j,k)**2)*ufluc_s(i,j,k)
          tke_s(i,j,k,1) = (tke_s(i,j,k,2)+tke_s(i,j,k,3)+tke_s(i,j,k,4))*0.5
          tke_s(i,j,k,5) = (ufluc_s(i,j,k)**2)*wfluc(i,j,k)
          tke_s(i,j,k,6) = (cfluc(i,j,k)**2)*ufluc_s(i,j,k)
          tke_s(i,j,k,7) = ufluc_s(i,j,k)*cfluc(i,j,k)*wfluc(i,j,k)
          tke_s(i,j,k,8) = (ufluc_s(i,j,k)**2)*cfluc(i,j,k)
           
          diff(i,j,k,2) = (ufluc_s(i,j,k)**2)
          diff(i,j,k,3) = (vfluc(i,j,k)**2)
          diff(i,j,k,4) = (wfluc(i,j,k)**2)
          diff(i,j,k,1) = (diff(i,j,k,2)+diff(i,j,k,3)+diff(i,j,k,4))*0.5
          diff(i,j,k,5) = (ufluc_s(i,j,k)*wfluc(i,j,k))
          diff(i,j,k,6) = (cfluc(i,j,k)**2)
          diff(i,j,k,7) = 0
          diff(i,j,k,8) = 0

          buo(i,j,k,:) = 0
          buo(i,j,k,1) = cfluc(i,j,k)*wfluc(i,j,k)
          buo(i,j,k,4) = cfluc(i,j,k)*wfluc(i,j,k)*2
          buo(i,j,k,5) = cfluc(i,j,k)*ufluc_s(i,j,k)
          buo(i,j,k,7) = cfluc(i,j,k)**2
          
          rad1(i,j,k,6) = qfluc(i,j,k)*cfluc(i,j,k)
          rad1(i,j,k,7) = qfluc(i,j,k)*wfluc(i,j,k)
          rad1(i,j,k,8) = qfluc(i,j,k)*ufluc_s(i,j,k)
          
          rd_s(i,j,k,6) = fxfluc(i,j,k)*cfluc(i,j,k)          
          rd_s(i,j,k,7) = fxfluc(i,j,k)*wfluc(i,j,k)          
          rd_s(i,j,k,8) = fxfluc(i,j,k)*ufluc_s(i,j,k)          
          
        end do
       end do
      end do

      call inter_s_c_m(imax,jmax*kmax/px*8,tke_s,tke_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px*8,pac_s,pac_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px*8,rd_s,rd_c,dr)

      call deriv_c_s_m(imax,jmax*kmax/px*8,tke_c,dtke_s,dr)
      call deriv_c_s_m(imax,jmax*kmax/px*8,pac_c,dpac_s,dr)
      call deriv_c_s_m(imax,jmax*kmax/px*8,rd_c,drd_s,dr)
      call deriv_s_c_m(imax,jmax*kmax/px*8,diff,diff_c,dr)

      do i=1,imax
       diff_c(i,:,:,:)=diff_c(i,:,:,:)*mr_c(i)
       hiums_c(i,:,:,:)=diff_c(i,:,:,:)
      enddo

      call deriv_c_s_m(imax,jmax*kmax/px*8,hiums_c,hiums_s,dr)

      do i=1,imax
       hiums_s(i,:,:,:)=hiums_s(i,:,:,:)*mr_s(i)
      enddo

      do rs=1,6
       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=1,imax
          difff(i,j,k,rs)=hiums_s(i,j,k,rs)
         enddo
        enddo
       enddo
      enddo
      difff2 = 0
      call produc(um,vm,wm,cm,prod,ufluc,vfluc,wfluc,cfluc)
      call dissp(eps,ps,epsrad,ufluc,vfluc,wfluc,pfluc,cfluc,difff2,fxfluc,fyfluc,fzfluc,dfyc_dt,dfzc_dz)

      do i=1,imax
       hiums_s(i,:,:,:)=difff2(i,:,:,:)
      enddo
      call deriv_s_c_m(imax,jmax*kmax/px*8,hiums_s,hiums_c,dr)

      do i=1,imax
       hiums_c(i,:,:,:)=hiums_c(i,:,:,:)*mr_c(i)
      enddo
      call inter_c_s_m(imax,jmax*kmax/px*8,hiums_c,hiums_s,dr)

      do rs=7,8
       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=1,imax
          difff(i,j,k,rs)=hiums_s(i,j,k,rs)
         enddo
        enddo
       enddo
      enddo

      buo=buo*Gr/(Re**2)
      do rs=1,5
       difff(:,:,:,rs)=difff(:,:,:,rs)/Re
       eps(:,:,:,rs)=eps(:,:,:,rs)/Re
      enddo
      do rs=6,8
       difff(:,:,:,rs)=difff(:,:,:,rs)/(Re*Pr)
       eps(:,:,:,rs)=eps(:,:,:,rs)/(Re*Pr)*(Pr+1.) 
       rad1(:,:,:,rs)=rad1(:,:,:,rs)/(Re*Pr*Pl)
      end do
      eps(:,:,:,6)=eps(:,:,:,6)/(Pr+1.)*2
      rad1(:,:,:,6)=rad1(:,:,:,6)*2
      drd_s=drd_s/(Re*Pr*Pl)
      drd_s(:,:,:,6) =drd_s(:,:,:,6)*2
      epsrad=epsrad/(Re*Pr*Pl)
      epsrad(:,:,:,6)=epsrad(:,:,:,6)*2

      do rs=1,8
       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=1,imax
       
          epsm(i,rs)=epsm(i,rs) - eps(i,j,k,rs)/(jmax*kmax*nf)
          diffm(i,rs)=diffm(i,rs) + difff(i,j,k,rs)/(jmax*kmax*nf)
          prodm(i,rs)=prodm(i,rs) - prod(i,j,k,rs)/(jmax*kmax*nf) 
          pdm(i,rs)=pdm(i,rs) - dpac_s(i,j,k,rs)*mr_s(i)/(jmax*kmax*nf)
          tdm(i,rs)=tdm(i,rs) - dtke_s(i,j,k,rs)*mr_s(i)/(jmax*kmax*nf)
          rdm(i,rs)=rdm(i,rs) - drd_s(i,j,k,rs)*mr_s(i)/(jmax*kmax*nf)
          psm(i,rs)=psm(i,rs) + ps(i,j,k,rs)/(jmax*kmax*nf)
          buom(i,rs)=buom(i,rs) + buo(i,j,k,rs)/(jmax*kmax*nf) 
          radm(i,rs)=radm(i,rs) - rad1(i,j,k,rs)/(jmax*kmax*nf)     
          rem(i,rs)=rem(i,rs) + epsrad(i,j,k,rs)/(jmax*kmax*nf)       
         enddo
        enddo 
       enddo
      enddo
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         dfyc(i) = dfyc(i)+dfyc_dt(i,j,k)/(jmax*kmax*nf)
         dfzc(i) = dfzc(i)+dfzc_dz(i,j,k)/(jmax*kmax*nf)
        enddo
       enddo
      enddo
      return
      end 

      subroutine produc(um,vm,wm,cm,prodd,utmp,vtmp,wtmp,ctmp)
      include 'par_post.txt'
      include 'common.txt'

      real utmp(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real ctmp(0:imax+1,jmax/p_row,kmax/p_col)
      real hiu_s(imax,jmax/p_row,kmax/p_col),hiu_c(0:imax,jmax/p_row,kmax/p_col)
      real hiv_s(imax,jmax/p_row,kmax/p_col),hiv_c(0:imax,jmax/p_row,kmax/p_col)
      real hiw_s(imax,jmax/p_row,kmax/p_col),hiw_c(0:imax,jmax/p_row,kmax/p_col)
      real hic_s(imax,jmax/p_row,kmax/p_col),hic_c(0:imax,jmax/p_row,kmax/p_col)
      real u_s(0:imax+1,jmax/p_row,kmax/p_col),u_c(0:imax+1,jmax/p_row,kmax/p_col)
      real v_s(0:imax+1,jmax/p_row,kmax/p_col),v_c(0:imax+1,jmax/p_row,kmax/p_col)
      real w_s(0:imax+1,jmax/p_row,kmax/p_col),w_c(0:imax+1,jmax/p_row,kmax/p_col)
      real c_s(0:imax+1,jmax/p_row,kmax/p_col),c_c(0:imax+1,jmax/p_row,kmax/p_col)

      real prodd(0:imax+1,jmax/p_row,kmax/p_col,8)
      real um(imax),vm(imax),wm(imax),wm_c(0:imax),dwm_s(imax),dwm_c(0:imax)
      real cm(imax),cm_c(0:imax),dcm_s(imax),dcm_c(0:imax)    
      integer rs


      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         hiv_s(i,j,k) = vtmp(i,j,k)
         hiw_s(i,j,k) = wtmp(i,j,k)
         hic_s(i,j,k) = ctmp(i,j,k)
        enddo
        do i=0,imax
         hiu_c(i,j,k) = utmp(i,j,k)
        enddo
       enddo
      enddo

      call inter_c_s_m(imax,jmax*kmax/px,hiu_c,hiu_s,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiv_s,hiv_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiw_s,hiw_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hic_s,hic_c,dr)
      call inter_s_c(imax,wm,wm_c,dr)
      call inter_s_c(imax,cm,cm_c,dr)
      call der1w_s_c6_m(imax,1,wm,dwm_c,dr)
      call der1c_s_c6_m(imax,1,cm,dcm_c,dr,Twb,Twt)
      dwm_c = dwm_c * mr_c
      dcm_c = dcm_c * mr_c
      call inter_c_s(imax,dwm_c,dwm_s,dr)
      call inter_c_s(imax,dcm_c,dcm_s,dr)
 
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         u_s(i,j,k) = hiu_s(i,j,k)
         v_s(i,j,k) = hiv_s(i,j,k)
         w_s(i,j,k) = hiw_s(i,j,k)
         c_s(i,j,k) = hic_s(i,j,k)
        enddo
        do i=0,imax
         u_c(i,j,k) = hiu_c(i,j,k)
         v_c(i,j,k) = hiv_c(i,j,k)
         w_c(i,j,k) = hiw_c(i,j,k)
         c_c(i,j,k) = hic_c(i,j,k)
        enddo
       enddo
      enddo
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
          prodd(i,j,k,1) = u_s(i,j,k)*w_s(i,j,k)*dwm_s(i)
          prodd(i,j,k,2) = 0.
          prodd(i,j,k,3) = 0.
          prodd(i,j,k,4) = 2.*u_s(i,j,k)*w_s(i,j,k)*dwm_s(i)
          prodd(i,j,k,5) = u_s(i,j,k)*u_s(i,j,k)*dwm_s(i)
          prodd(i,j,k,6) = 2.*u_s(i,j,k)*c_s(i,j,k)*dcm_s(i)
          prodd(i,j,k,7) = u_s(i,j,k)*c_s(i,j,k)*dwm_s(i)+u_s(i,j,k)*w_s(i,j,k)*dcm_s(i)
          prodd(i,j,k,8) = u_s(i,j,k)*u_s(i,j,k)*dcm_s(i)
        enddo
       enddo
      enddo
      end
 
 
!------------------------------------- LOAD SUBROUTINE

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
      real rho_c(0:imax,jmax/p_row,kmax/p_col)
      real rho_s(0:i1,jmax/p_row,kmax/p_col)
      real pi
      integer istap
      character*5 cha
      character*5 cha2
      call cnvstr(nnrank,cha)
      call cnvstr(istap,cha2)
      ptmp = 0
      pi = 4.*atan(1.0)

      if (ini.eq.0) then
      if (nrank.eq.0)  write (*,*) cha2,imax,kmax,jmax
       call decomp_2d_read_one(1,unew,'DATA/u.'//cha2//'.dns')
       call decomp_2d_read_one(1,vnew,'DATA/v.'//cha2//'.dns')
       call decomp_2d_read_one(1,wnew,'DATA/w.'//cha2//'.dns')
       call decomp_2d_read_one(1,cnew,'DATA/c.'//cha2//'.dns')
       call decomp_2d_read_one(1,ptmp,'DATA/p.'//cha2//'.dns')
       call decomp_2d_read_one(1,qrad,'DATA/q.'//cha2//'.dns')
       call decomp_2d_read_one(1,Gnew,'DATA/G.'//cha2//'.dns')
       call decomp_2d_read_one(1,knew,'DATA/k.'//cha2//'.dns')
       call decomp_2d_read_one(1,vcor,'DATA/vc.'//cha2//'.dns')
       call decomp_2d_read_one(1,fxnew,'DATA/k.'//cha2//'.dns')
       call decomp_2d_read_one(1,fynew,'DATA/k.'//cha2//'.dns')
       call decomp_2d_read_one(1,fznew,'DATA/k.'//cha2//'.dns')
      endif
       enew=4 * (cnew/Tplus+1)**4
       Gnew = Gnew/knew
       fxnew=fxnew/pi
       fynew=fynew/pi
       fznew=fznew/pi
       do i=1,imax
        p(i,:,:) = ptmp(i,:,:)
       enddo
       call stateIG() 
       do i=1,imax
        do j=1,jmax/p_row
         do k=1,kmax/p_col
          renew(i,j,k) = cnew(i,j,k)*rhonew(i,j,k)
          phitnew(i,j,k) = rhonew(i,j,k)*vnew(i,j,k)
          phiznew(i,j,k) = rhonew(i,j,k)*wnew(i,j,k)
         enddo
        enddo
       enddo
       rho_s=rhonew
       call inter_s_c_m(imax,jmax*kmax/px,rho_s,rho_c,dr)
       do i=1,imax
        do j=1,jmax/p_row
         do k=1,kmax/p_col
          phirnew(i,j,k) = rho_c(i,j,k)*unew(i,j,k)
         enddo
        enddo
       enddo
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
      end



      subroutine parsum(aam,bbm,ccm,ddm,eem,ffm,ggm)
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real mpi_sm(imax)
      real aam(imax),bbm(imax),ccm(imax),ddm(imax),eem(imax),ffm(imax),ggm(imax)

      call mpi_allreduce(aam,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      aam = mpi_sm

      call mpi_allreduce(bbm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      bbm = mpi_sm

      call mpi_allreduce(ccm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ccm = mpi_sm

      call mpi_allreduce(ddm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ddm = mpi_sm

      call mpi_allreduce(eem,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      eem = mpi_sm

      call mpi_allreduce(ffm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ffm = mpi_sm

      call mpi_allreduce(ggm,mpi_sm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ggm = mpi_sm

      end


      subroutine parsum5(aam,bbm,ccm,ddm,eem,ffm,ggm)
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real mpi_sm(imax,8)
      real aam(imax,8),bbm(imax,8),ccm(imax,8),ddm(imax,8),eem(imax,8),ffm(imax,8),ggm(imax,8)

      call mpi_allreduce(aam,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      aam = mpi_sm

      call mpi_allreduce(bbm,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      bbm = mpi_sm

      call mpi_allreduce(ccm,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ccm = mpi_sm

      call mpi_allreduce(ddm,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ddm = mpi_sm

      call mpi_allreduce(eem,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      eem = mpi_sm

      call mpi_allreduce(ffm,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ffm = mpi_sm

      call mpi_allreduce(ggm,mpi_sm,imax*8,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ggm = mpi_sm

      end

      subroutine parsum4(aam,bbm)
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real mpi_sm(imax,4)
      real aam(imax,4),bbm(imax,4)
      call mpi_allreduce(aam,mpi_sm,imax*4,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      aam = mpi_sm
      call mpi_allreduce(bbm,mpi_sm,imax*4,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      bbm = mpi_sm
      end

      subroutine parsum15(aam)
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real mpi_sm(imax,15)
      real aam(imax,15)
      call mpi_allreduce(aam,mpi_sm,imax*15,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      aam = mpi_sm
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


      subroutine dissp(epss,ps,epsrad,utmp,vtmp,wtmp,ptmp,ctmp,difff,fxtmp,fytmp,fztmp,dfyc_dt,dfzc_dz)
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'

      real d11,d12,d13,d21,d22,d23,d31,d32,d33,dc1,dc2,dc3
      real s11,s12,s13,s21,s22,s23,s31,s32,s33
      real utmp(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real ctmp(0:imax+1,jmax/p_row,kmax/p_col)
      real ptmp(0:imax+1,jmax/p_row,kmax/p_col)
      real fxtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real fytmp(0:imax+1,jmax/p_row,kmax/p_col)
      real fztmp(0:imax+1,jmax/p_row,kmax/p_col)
      real himu_s(imax,jmax/p_row,kmax/p_col),himu_c(0:imax,jmax/p_row,kmax/p_col)
      real himum_s(imax,jmax/p_row,kmax/p_col),himum_c(0:imax,jmax/p_row,kmax/p_col)
      real hiu_s(imax,jmax/p_row,kmax/p_col),hiu_c(0:imax,jmax/p_row,kmax/p_col)
      real hiv_s(imax,jmax/p_row,kmax/p_col),hiv_c(0:imax,jmax/p_row,kmax/p_col)
      real hiw_s(imax,jmax/p_row,kmax/p_col),hiw_c(0:imax,jmax/p_row,kmax/p_col)
      real hic_s(imax,jmax/p_row,kmax/p_col),hic_c(0:imax,jmax/p_row,kmax/p_col)
      real u_s(0:imax+1,jmax/p_row,kmax/p_col),u_c(0:imax+1,jmax/p_row,kmax/p_col)
      real v_s(0:imax+1,jmax/p_row,kmax/p_col),v_c(0:imax+1,jmax/p_row,kmax/p_col)
      real w_s(0:imax+1,jmax/p_row,kmax/p_col),w_c(0:imax+1,jmax/p_row,kmax/p_col)
      real c_s(0:imax+1,jmax/p_row,kmax/p_col),c_c(0:imax+1,jmax/p_row,kmax/p_col)
      real fx_s(0:imax+1,jmax/p_row,kmax/p_col)
      real fy_s(0:imax+1,jmax/p_row,kmax/p_col)
      real fz_s(0:imax+1,jmax/p_row,kmax/p_col)
      real fzc_s(0:imax+1,jmax/p_row,kmax/p_col)
      real fyc_s(0:imax+1,jmax/p_row,kmax/p_col)
      real dhiu_s(imax,jmax/p_row,kmax/p_col),dhiu_c(0:imax,jmax/p_row,kmax/p_col)
      real dhiv_s(imax,jmax/p_row,kmax/p_col),dhiv_c(0:imax,jmax/p_row,kmax/p_col)
      real dhiw_s(imax,jmax/p_row,kmax/p_col),dhiw_c(0:imax,jmax/p_row,kmax/p_col)
      real dhic_s(imax,jmax/p_row,kmax/p_col),dhic_c(0:imax,jmax/p_row,kmax/p_col)
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax)
      real tmp (0:imx,jmax,kmax/p_col),tmpf (0:imx,jmax,kmax/p_col)
      real tmp2(0:imx,jmax,kmax/p_col),tmp2f(0:imx,jmax,kmax/p_col)
      real tmp3(0:imx,jmax/p_col,kmax),tmp3f(0:imx,jmax/p_col,kmax)
      real du_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real du_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dfyc_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dfzc_dz(0:imax+1,jmax/p_row,kmax/p_col)

      integer rs

      real epss(0:imax+1,jmax/p_row,kmax/p_col,8)
      real epsrad(0:imax+1,jmax/p_row,kmax/p_col,8)
      real difff(0:imax+1,jmax/p_row,kmax/p_col,8)
      real ps(0:imax+1,jmax/p_row,kmax/p_col,8)

      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         hiv_s(i,j,k) = vtmp(i,j,k)
         hiw_s(i,j,k) = wtmp(i,j,k)
         hic_s(i,j,k) = ctmp(i,j,k)
         fx_s(i,j,k) = fxtmp(i,j,k)
         fy_s(i,j,k) = fytmp(i,j,k)
         fz_s(i,j,k) = fztmp(i,j,k)
         fyc_s(i,j,k) = fytmp(i,j,k)*ctmp(i,j,k)
         fzc_s(i,j,k) = fztmp(i,j,k)*ctmp(i,j,k)
        enddo
        do i=0,imax
         hiu_c(i,j,k) = utmp(i,j,k)
        enddo
       enddo
      enddo


      call inter_c_s_m(imax,jmax*kmax/px,hiu_c,hiu_s,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiv_s,hiv_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiw_s,hiw_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hic_s,hic_c,dr)

      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         u_s(i,j,k) = hiu_s(i,j,k)
         v_s(i,j,k) = hiv_s(i,j,k)
         w_s(i,j,k) = hiw_s(i,j,k)
         c_s(i,j,k) = hic_s(i,j,k)
        enddo
        do i=0,imax
         u_c(i,j,k)  = hiu_c(i,j,k)
         v_c(i,j,k)  = hiv_c(i,j,k)
         w_c(i,j,k)  = hiw_c(i,j,k)
         c_c(i,j,k)  = hic_c(i,j,k)
        enddo
       enddo
      enddo


       call deriv_c_s_m(imax,jmax*kmax/px,hiu_c,dhiu_s,dr)
       call deriv_c_s_m(imax,jmax*kmax/px,hiv_c,dhiv_s,dr)
       call deriv_c_s_m(imax,jmax*kmax/px,hiw_c,dhiw_s,dr)
       call deriv_c_s_m(imax,jmax*kmax/px,hic_c,dhic_s,dr)

! -------------------------------------------------------------------

       call transpose_x_to_y(u_s,tmp)

       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft (j)=tmp (i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp (i,j,k)=dft (j)
          enddo
         enddo
      enddo
      call transpose_y_to_x(tmp,du_dt)
!-------------------------------------------------------------------
      call transpose_x_to_y(u_s,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
         do j=1,jmax/p_col
           do k=1,kmax
             fz (k)=tmp3 (i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3 (i,j,k)=dfz(k)
           enddo
             enddo
            enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,du_dz)
!------------------------------------------------------------------
!-------------------------------------------------------------------

        call transpose_x_to_y(v_s,tmp)
       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft(j)=tmp(i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp (i,j,k)=dft (j)
          enddo
         enddo
      enddo

      call transpose_y_to_x(tmp,dv_dt)
!---------------------------------------------------------------------

      call transpose_x_to_y(v_s,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
            do j=1,jmax/p_col
           do k=1,kmax
             fz (k)=tmp3 (i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3 (i,j,k)=dfz(k)
           enddo
             enddo
            enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,dv_dz)

!----------------------------------------------------------------------
!----------------------------------------------------------------------

       call transpose_x_to_y(w_s,tmp)

       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft(j)=tmp(i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp (i,j,k)=dft (j)
          enddo
         enddo
      enddo

      call transpose_y_to_x(tmp,dw_dt)
!-------------------------------------------------------------------------

      call transpose_x_to_y(w_s,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
            do j=1,jmax/p_col
           do k=1,kmax
             fz (k)=tmp3 (i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3 (i,j,k)=dfz(k)
           enddo
             enddo
            enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,dw_dz)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

       call transpose_x_to_y(c_s,tmp)

       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft(j)=tmp(i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp(i,j,k)=dft(j)
          enddo
         enddo
      enddo

      call transpose_y_to_x(tmp,dc_dt)
!-------------------------------------------------------------------------
 
      call transpose_x_to_y(c_s,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
            do j=1,jmax/p_col
           do k=1,kmax
             fz(k)=tmp3(i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3(i,j,k)=dfz(k)
           enddo
             enddo
            enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,dc_dz)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

       call transpose_x_to_y(fyc_s,tmp)

       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft(j)=tmp(i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp(i,j,k)=dft(j)
          enddo
         enddo
      enddo

      call transpose_y_to_x(tmp,dfyc_dt)
!-------------------------------------------------------------------------
 
      call transpose_x_to_y(fzc_s,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
            do j=1,jmax/p_col
           do k=1,kmax
             fz(k)=tmp3(i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3(i,j,k)=dfz(k)
           enddo
             enddo
            enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,dfzc_dz)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax

         d11= dhiu_s(i,j,k)*mr_s(i)
         d12= du_dt(i,j,k)
         d13= du_dz(i,j,k)

         d21= dhiv_s(i,j,k)*mr_s(i)
         d22= dv_dt(i,j,k)
         d23= dv_dz(i,j,k)

         d31= dhiw_s(i,j,k)*mr_s(i)
         d32= dw_dt(i,j,k)
         d33= dw_dz(i,j,k)

         dc1= dhic_s(i,j,k)*mr_s(i)
         dc2= dc_dt(i,j,k)
         dc3= dc_dz(i,j,k)
        
        epss(i,j,k,1)= d11**2+d12**2+d13**2+
     +                 d21**2+d22**2+d23**2+
     +                 d31**2+d32**2+d33**2
        ps(i,j,k,1)= 0

        epss(i,j,k,2)= 2*(d11**2+d12**2+d13**2) 

        ps(i,j,k,2)= 2.*ptmp(i,j,k)*(d11)

        epss(i,j,k,3)= 2*(d21**2+d22**2+d23**2) 

        ps(i,j,k,3)= 2.*ptmp(i,j,k)*(d22)

        epss(i,j,k,4)= 2*(d31**2+d32**2+d33**2) 

        ps(i,j,k,4)= 2.*ptmp(i,j,k)*(d33)

        epss(i,j,k,5)= 2*(d11*d31+d12*d32+d13*d33) 

        ps(i,j,k,5)= ptmp(i,j,k)*(d13+d31)

        epss(i,j,k,6)= dc1**2+dc2**2+dc3**2

        ps(i,j,k,6)= 0

        epss(i,j,k,7)= dc1*d31+dc2*d32+dc3*d33

        ps(i,j,k,7)= ptmp(i,j,k)*(dc3)

        epss(i,j,k,8)= dc1*d11+dc2*d12+dc3*d13

        ps(i,j,k,8)= ptmp(i,j,k)*(dc1)

        difff(i,j,k,7)= Pr*c_s(i,j,k)*d31+w_s(i,j,k)*dc1
 
        difff(i,j,k,8)= Pr*c_s(i,j,k)*d11+u_s(i,j,k)*dc1

        epsrad(i,j,k,6)= fx_s(i,j,k)*dc1+fy_s(i,j,k)*dc2+fz_s(i,j,k)*dc3

        epsrad(i,j,k,7)= fx_s(i,j,k)*d31+fy_s(i,j,k)*d32+fz_s(i,j,k)*d33

        epsrad(i,j,k,8)= fx_s(i,j,k)*d11+fy_s(i,j,k)*d12+fz_s(i,j,k)*d13

        enddo
       enddo
      enddo
      
      end
     
      subroutine spectrazcalc(dcc,ufluc,ufluc_s,vfluc,wfluc,cfluc,pfluc,qfluc,efluc,Pefluc,Pafluc,kfluc,fxfluc,fyfluc,fzfluc,
     +           Gfluc,speczuu,speczGG,
     +           speczvv,speczww,speczuw,speczuc,speczwc,speczqq,speczqc,speczcc,speczec,speczee,speczkk,speczfxc,specbudz,
     +               nfiles)

      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real ufluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc_s(0:imax+1,jmax/p_row, kmax/p_col)
      real vfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real pfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real cfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real efluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Gfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pefluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pafluc(0:imax+1,jmax/p_row, kmax/p_col)
      real qfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real kfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real fxfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fyfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fzfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ucfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real vcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dvcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dwcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ducfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real specuu1(0:imx,kmax),specvv1(0:imx,kmax),specww1(0:imx,kmax)
      real specuw1(0:imx,kmax),specpp1(0:imx,kmax),specqq1(0:imx,kmax)
      real speccc1(0:imx,kmax),specuc1(0:imx,kmax),specwc1(0:imx,kmax)
      real specqc1(0:imx,kmax),specec1(0:imx,kmax),speckk1(0:imx,kmax)
      real specee1(0:imx,kmax)
      real specGG1(0:imx,kmax)
      real specfxc1(0:imx,kmax)
      real speczuu(0:imx,0:kmax/2),speczvv(0:imx,0:kmax/2),speczww(0:imx,0:kmax/2)
      real speczuw(0:imx,0:kmax/2),speczpp(0:imx,0:kmax/2),speczqq(0:imx,0:kmax/2)
      real speczcc(0:imx,0:kmax/2),speczuc(0:imx,0:kmax/2),speczwc(0:imx,0:kmax/2)
      real speczqc(0:imx,0:kmax/2),speczec(0:imx,0:kmax/2),speczkk(0:imx,0:kmax/2)
      real speczee(0:imx,0:kmax/2)
      real speczGG(0:imx,0:kmax/2)
      real speczfxc(0:imx,0:kmax/2)
      real specep1(0:imx,kmax),specep2(0:imx,kmax),spectp1(0:imx,kmax)
      real speczep1(0:imx,0:kmax/2),speczep2(0:imx,0:kmax/2),specztp(0:imx,0:kmax/2)
      real specty1(0:imx,kmax)
      real speczty(0:imx,0:kmax/2)
      real specsw1(0:imx,kmax),specer1(0:imx,kmax),specer2(0:imx,kmax)
      real speczsw(0:imx,0:kmax/2),speczer1(0:imx,0:kmax/2),speczer2(0:imx,0:kmax/2)
      real specfzc1(0:imx,kmax)
      real speczfzc(0:imx,0:kmax/2)
      real specbudz(0:imx,0:kmax/2,14)
      real specPe1(0:imx,kmax),specPa1(0:imx,kmax)
      real speczPe(0:imx,0:kmax/2),speczPa(0:imx,0:kmax/2)

      real interp_c(0:imax,jmax/p_row,kmax/p_col)
      real interp_s(imax,jmax/p_row,kmax/p_col)
      real deriv_s(imax,jmax/p_row,kmax/p_col)
      real zint_c(0:imax,0:kmax/2)
      real zint_s(imax,0:kmax/2)
      real zder_s(imax,0:kmax/2)
      real zder2_s(imax,0:kmax/2)

      real dcc(imax)
      integer nfiles

      ducfluc=0
      dcfluc=0
      ucfluc=ufluc_s*cfluc
      wcfluc=wfluc*cfluc
      vcfluc=vfluc*cfluc
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dcfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      do i=1,imax
       interp_s(i,:,:)=ucfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       ducfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      do i=1,imax
       interp_s(i,:,:)=vcfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dvcfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      call yzderiv(cfluc,dc_dt,dc_dz)

      call spectraz(ufluc,ufluc,specuu1,nfiles)
      call spectraz(vfluc,vfluc,specvv1,nfiles)
      call spectraz(wfluc,wfluc,specww1,nfiles)
      call spectraz(cfluc,cfluc,speccc1,nfiles)
      call spectraz(efluc,cfluc,specec1,nfiles)
      call spectraz(pfluc,pfluc,specpp1,nfiles)
      call spectraz(qfluc,qfluc,specqq1,nfiles)
      call spectraz(ufluc_s,wfluc,specuw1,nfiles)
      call spectraz(ufluc_s,cfluc,specuc1,nfiles)
      call spectraz(wfluc,cfluc,specwc1,nfiles)
      call spectraz(qfluc,cfluc,specqc1,nfiles)
      call spectraz(kfluc,kfluc,speckk1,nfiles)
      call spectraz(efluc,efluc,specee1,nfiles)
      call spectraz(Gfluc,Gfluc,specGG1,nfiles)
      call spectraz(Pefluc,cfluc,specPe1,nfiles)
      call spectraz(Pafluc,cfluc,specPa1,nfiles)
      call spectraz(fxfluc,cfluc,specfxc1,nfiles)
      call spectrazi(fzfluc,cfluc,specfzc1,nfiles)
      call spectraz(dc_dt,dc_dt,specep1,nfiles)
      call spectraz(dcfluc,dcfluc,specep2,nfiles)
      call spectraz(ducfluc,cfluc,spectp1,nfiles)
      call spectraz(dvcfluc,cfluc,specty1,nfiles)
      call spectrazi(cfluc,wcfluc,specsw1,nfiles)
      call spectraz(fxfluc,dcfluc,specer1,nfiles)
      call spectraz(fyfluc,dc_dt,specer2,nfiles)


      do i=0,imx
         do k=2,kmax-2,2
        speczuu(i,k/2)=speczuu(i,k/2)+(specuu1(i,k)+specuu1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczvv(i,k/2)=speczvv(i,k/2)+(specvv1(i,k)+specvv1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczww(i,k/2)=speczww(i,k/2)+(specww1(i,k)+specww1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczuw(i,k/2)=speczuw(i,k/2)+(specuw1(i,k)+specuw1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczpp(i,k/2)=speczpp(i,k/2)+(specpp1(i,k)+specpp1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczqq(i,k/2)=speczqq(i,k/2)+(specqq1(i,k)+specqq1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczcc(i,k/2)=speczcc(i,k/2)+(speccc1(i,k)+speccc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczec(i,k/2)=speczec(i,k/2)+(specec1(i,k)+specec1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczee(i,k/2)=speczee(i,k/2)+(specee1(i,k)+specee1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczGG(i,k/2)=speczGG(i,k/2)+(specGG1(i,k)+specGG1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczPe(i,k/2)=speczPe(i,k/2)+(specPe1(i,k)+specPe1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczPa(i,k/2)=speczPa(i,k/2)+(specPa1(i,k)+specPa1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczuc(i,k/2)=speczuc(i,k/2)+(specuc1(i,k)+specuc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczwc(i,k/2)=speczwc(i,k/2)+(specwc1(i,k)+specwc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczqc(i,k/2)=speczqc(i,k/2)+(specqc1(i,k)+specqc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczkk(i,k/2)=speczkk(i,k/2)+(speckk1(i,k)+speckk1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczfxc(i,k/2)=speczfxc(i,k/2)+(specfxc1(i,k)+specfxc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczfzc(i,k/2)=speczfzc(i,k/2)+(specfzc1(i,k)+specfzc1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczsw(i,k/2)=speczsw(i,k/2)+(specsw1(i,k)+specsw1(i,k+1))/(jmax*nfiles)*2./(kmax)
        specztp(i,k/2)=specztp(i,k/2)+(spectp1(i,k)+spectp1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczty(i,k/2)=speczty(i,k/2)+(specty1(i,k)+specty1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczep1(i,k/2)=speczep1(i,k/2)+(specep1(i,k)+specep1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczep2(i,k/2)=speczep2(i,k/2)+(specep2(i,k)+specep2(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczer1(i,k/2)=speczer1(i,k/2)+(specer1(i,k)+specer1(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczer2(i,k/2)=speczer2(i,k/2)+(specer2(i,k)+specer2(i,k+1))/(jmax*nfiles)*2./(kmax)
         enddo

        speczuu(i,0)=speczuu(i,0)+(specuu1(i,1))/(jmax*kmax*nfiles)
        speczvv(i,0)=speczvv(i,0)+(specvv1(i,1))/(jmax*kmax*nfiles)
        speczww(i,0)=speczww(i,0)+(specww1(i,1))/(jmax*kmax*nfiles)
        speczuw(i,0)=speczuw(i,0)+(specuw1(i,1))/(jmax*kmax*nfiles)
        speczpp(i,0)=speczpp(i,0)+(specpp1(i,1))/(jmax*kmax*nfiles)
        speczqq(i,0)=speczqq(i,0)+(specqq1(i,1))/(jmax*nfiles*kmax)
        speczcc(i,0)=speczcc(i,0)+(speccc1(i,1))/(jmax*nfiles*kmax)
        speczec(i,0)=speczec(i,0)+(specec1(i,1))/(jmax*nfiles*kmax)
        speczee(i,0)=speczee(i,0)+(specee1(i,1))/(jmax*nfiles*kmax)
        speczGG(i,0)=speczGG(i,0)+(specGG1(i,1))/(jmax*nfiles*kmax)
        speczPe(i,0)=speczPe(i,0)+(specPe1(i,1))/(jmax*nfiles*kmax)
        speczPa(i,0)=speczPa(i,0)+(specPa1(i,1))/(jmax*nfiles*kmax)
        speczuc(i,0)=speczuc(i,0)+(specuc1(i,1))/(jmax*nfiles*kmax)
        speczwc(i,0)=speczwc(i,0)+(specwc1(i,1))/(jmax*nfiles*kmax)
        speczqc(i,0)=speczqc(i,0)+(specqc1(i,1))/(jmax*nfiles*kmax)
        speczkk(i,0)=speczkk(i,0)+(speckk1(i,1))/(jmax*nfiles*kmax)
        speczfxc(i,0)=speczfxc(i,0)+(specfxc1(i,1))/(jmax*nfiles*kmax)
        speczfzc(i,0)=speczfzc(i,0)+(specfzc1(i,1))/(jmax*nfiles*kmax)
        speczsw(i,0)=speczsw(i,0)+(specsw1(i,1))/(jmax*nfiles*kmax)
        specztp(i,0)=specztp(i,0)+(spectp1(i,1))/(jmax*nfiles*kmax)
        speczty(i,0)=speczty(i,0)+(specty1(i,1))/(jmax*nfiles*kmax)
        speczep1(i,0)=speczep1(i,0)+(specep1(i,1))/(jmax*nfiles*kmax)
        speczep2(i,0)=speczep2(i,0)+(specep2(i,1))/(jmax*nfiles*kmax)
        speczer1(i,0)=speczer1(i,0)+(specer1(i,1))/(jmax*nfiles*kmax)
        speczer2(i,0)=speczer2(i,0)+(specer2(i,1))/(jmax*nfiles*kmax)

        speczuu(i,kmax/2)=speczuu(i,kmax/2)+(specuu1(i,kmax))/(jmax*kmax*nfiles)*2.
        speczvv(i,kmax/2)=speczvv(i,kmax/2)+(specvv1(i,kmax))/(jmax*kmax*nfiles)*2.
        speczww(i,kmax/2)=speczww(i,kmax/2)+(specww1(i,kmax))/(jmax*kmax*nfiles)*2.
        speczuw(i,kmax/2)=speczuw(i,kmax/2)+(specuw1(i,kmax))/(jmax*kmax*nfiles)*2.
        speczpp(i,kmax/2)=speczpp(i,kmax/2)+(specpp1(i,kmax))/(jmax*kmax*nfiles)*2.
        speczqq(i,kmax/2)=speczqq(i,kmax/2)+(specqq1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczcc(i,kmax/2)=speczcc(i,kmax/2)+(speccc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczec(i,kmax/2)=speczec(i,kmax/2)+(specec1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczee(i,kmax/2)=speczee(i,kmax/2)+(specee1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczGG(i,kmax/2)=speczGG(i,kmax/2)+(specGG1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczPe(i,kmax/2)=speczPe(i,kmax/2)+(specPe1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczPa(i,kmax/2)=speczPa(i,kmax/2)+(specPa1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczuc(i,kmax/2)=speczuc(i,kmax/2)+(specuc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczwc(i,kmax/2)=speczwc(i,kmax/2)+(specwc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczqc(i,kmax/2)=speczqc(i,kmax/2)+(specqc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczkk(i,kmax/2)=speczkk(i,kmax/2)+(speckk1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczfxc(i,kmax/2)=speczfxc(i,kmax/2)+(specfxc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczfzc(i,kmax/2)=speczfzc(i,kmax/2)+(specfzc1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczsw(i,kmax/2)=speczsw(i,kmax/2)+(specsw1(i,kmax))/(jmax*nfiles*kmax)*2.
        specztp(i,kmax/2)=specztp(i,kmax/2)+(spectp1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczty(i,kmax/2)=speczty(i,kmax/2)+(specty1(i,kmax))/(jmax*nfiles*kmax)*2.
        speczep1(i,kmax/2)=speczep1(i,kmax/2)+(specep1(i,kmax))/(jmax*nfiles)*2./(kmax)
        speczep2(i,kmax/2)=speczep2(i,kmax/2)+(specep2(i,kmax))/(jmax*nfiles)*2./(kmax)
        speczer1(i,kmax/2)=speczer1(i,kmax/2)+(specer1(i,kmax))/(jmax*nfiles)*2./(kmax)
        speczer2(i,kmax/2)=speczer2(i,kmax/2)+(specer2(i,kmax))/(jmax*nfiles)*2./(kmax)

       enddo

      do i=1,imax
       zint_s(i,:)=speczcc(i,:)
      enddo
      call deriv_s_c_m(imax,kmax/2+1,zint_s,zint_c,dr)
      do k=0,kmax/2
       zint_c(:,k)=zint_c(:,k)*mr_c(:)
      enddo
      call deriv_c_s_m(imax,kmax/2+1,zint_c,zder_s,dr)
      do k=0,kmax/2
       zder_s(:,k)=zder_s(:,k)*mr_s(:)
      enddo
      zint_s=0.
      zint_c=0.
      zder2_s=0.
      do i=1,imax
       zint_s(i,:)=speczfxc(i,:)
      enddo
      call deriv_s_c_m(imax,kmax/2+1,zint_s,zint_c,dr)
       do k=0,kmax/2
        zint_c(:,k)=zint_c(:,k)*mr_c(:)
       enddo
      call inter_c_s_m(imax,kmax/2+1,zint_c,zder2_s,dr)

      do i=1,imx
       do k=0,kmax/2
        specbudz(i,k,1)=-2*speczuc(i,k)*dcc(i)
       enddo
       do k=0,kmax/2
        specbudz(i,k,2)=(speczsw(i,k)*(k*2./4.))
       enddo
       specbudz(i,:,3)=specztp(i,:)
       specbudz(i,:,4)=speczty(i,:)
       do k=0,kmax/2
        specbudz(i,k,5)=-2*speczcc(i,k)*(((k-1)*2./4.)**2.)/(Re*Pr)
       enddo
       specbudz(i,:,6)=-(speczep1(i,:))*2./(Re*Pr)
       specbudz(i,:,7)=-(speczep2(i,:))*2./(Re*Pr)
       specbudz(i,:,8)=(zder_s(i,:))/(Re*Pr)
       specbudz(i,:,9)=-(zder2_s(i,:))*2/(Re*Pr*Pl)
       specbudz(i,:,10)=(speczer1(i,:))*2/(Re*Pr*Pl)
       specbudz(i,:,11)=(speczer2(i,:))*2/(Re*Pr*Pl)
       do k=0,kmax/2
        specbudz(i,k,12)=(speczfzc(i,k)*(k*2./4.))*2/(Re*Pr*Pl)
       enddo
       specbudz(i,:,13)=-(speczPe(i,:))*2/(Re*Pr*Pl)
       specbudz(i,:,14)=(speczPa(i,:))*2/(Re*Pr*Pl)
      enddo

      return
      end
 
      subroutine spectraycalc(dcc,ufluc,ufluc_s,vfluc,wfluc,cfluc,pfluc,qfluc,efluc,Pefluc,Pafluc,kfluc,fxfluc,fyfluc,fzfluc,
     +             Gfluc,specyuu,specyGG,
     +             specyvv,specyww,specyuw,specyuc,specywc,specyqq,specyqc,specycc,specyec,specyee,specykk,specyfxc,specbudy,
     +               nfiles)

      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real ufluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ufluc_s(0:imax+1,jmax/p_row, kmax/p_col)
      real vfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real pfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real cfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real efluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Gfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pefluc(0:imax+1,jmax/p_row, kmax/p_col)
      real Pafluc(0:imax+1,jmax/p_row, kmax/p_col)
      real qfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real kfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real fxfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fyfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real fzfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ucfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real wcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real vcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dvcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dwcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real ducfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real specuu1(0:imx,jmax),specvv1(0:imx,jmax),specww1(0:imx,jmax)
      real specuw1(0:imx,jmax),specpp1(0:imx,jmax),specqq1(0:imx,jmax)
      real speccc1(0:imx,jmax),specuc1(0:imx,jmax),specwc1(0:imx,jmax)
      real specqc1(0:imx,jmax),specec1(0:imx,jmax),speckk1(0:imx,jmax)
      real specee1(0:imx,jmax)
      real specGG1(0:imx,jmax)
      real specfxc1(0:imx,jmax)
      real specyuu(0:imx,0:jmax/2),specyvv(0:imx,0:jmax/2),specyww(0:imx,0:jmax/2)
      real specyuw(0:imx,0:jmax/2),specypp(0:imx,0:jmax/2),specyqq(0:imx,0:jmax/2)
      real specycc(0:imx,0:jmax/2),specyuc(0:imx,0:jmax/2),specywc(0:imx,0:jmax/2)
      real specyqc(0:imx,0:jmax/2),specyec(0:imx,0:jmax/2),specykk(0:imx,0:jmax/2)
      real specyee(0:imx,0:jmax/2)
      real specyGG(0:imx,0:jmax/2)
      real specyfxc(0:imx,0:jmax/2)
      real specep1(0:imx,jmax),specep2(0:imx,jmax),spectp1(0:imx,jmax)
      real specyep1(0:imx,0:jmax/2),specyep2(0:imx,0:jmax/2),specytp(0:imx,0:jmax/2)
      real spectz1(0:imx,jmax)
      real specytz(0:imx,0:jmax/2)
      real specsv1(0:imx,jmax),specer1(0:imx,jmax),specer2(0:imx,jmax)
      real specysv(0:imx,0:jmax/2),specyer1(0:imx,0:jmax/2),specyer2(0:imx,0:jmax/2)
      real specfyc1(0:imx,jmax)
      real specyfyc(0:imx,0:jmax/2)
      real specbudy(0:imx,0:jmax/2,14),dcc(imax)
      real specPe1(0:imx,jmax),specPa1(0:imx,jmax)
      real specyPe(0:imx,0:jmax/2),specyPa(0:imx,0:jmax/2)
      integer nfiles

      real interp_c(0:imax,jmax/p_row,kmax/p_col)
      real interp_s(imax,jmax/p_row,kmax/p_col)
      real deriv_s(imax,jmax/p_row,kmax/p_col)
      real yint_c(0:imax,0:jmax/2)
      real yint_s(imax,0:jmax/2)
      real yder_s(imax,0:jmax/2)
      real yder2_s(imax,0:jmax/2)

      ducfluc=0
      dcfluc=0
      ucfluc=ufluc_s*cfluc
      wcfluc=wfluc*cfluc
      vcfluc=vfluc*cfluc
      do i=1,imax
       interp_s(i,:,:)=cfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dcfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      do i=1,imax
       interp_s(i,:,:)=ucfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       ducfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      do i=1,imax
       interp_s(i,:,:)=wcfluc(i,:,:)
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        interp_c(:,j,k)=interp_c(:,j,k)*mr_c(:)
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,interp_c,deriv_s,dr)
      do i=1,imax
       dwcfluc(i,:,:)=deriv_s(i,:,:)
      enddo

      call yzderiv(cfluc,dc_dt,dc_dz)

      call spectray(ufluc,ufluc,specuu1,nfiles)
      call spectray(vfluc,vfluc,specvv1,nfiles)
      call spectray(wfluc,wfluc,specww1,nfiles)
      call spectray(cfluc,cfluc,speccc1,nfiles)
      call spectray(efluc,cfluc,specec1,nfiles)
      call spectray(pfluc,pfluc,specpp1,nfiles)
      call spectray(qfluc,qfluc,specqq1,nfiles)
      call spectray(ufluc_s,wfluc,specuw1,nfiles)
      call spectray(ufluc_s,cfluc,specuc1,nfiles)
      call spectray(wfluc,cfluc,specwc1,nfiles)
      call spectray(qfluc,cfluc,specqc1,nfiles)
      call spectray(kfluc,kfluc,speckk1,nfiles)
      call spectray(efluc,efluc,specee1,nfiles)
      call spectray(Gfluc,Gfluc,specGG1,nfiles)
      call spectray(Pefluc,cfluc,specPe1,nfiles)
      call spectray(Pafluc,cfluc,specPa1,nfiles)
      call spectray(fxfluc,cfluc,specfxc1,nfiles)
      call spectrayi(fyfluc,cfluc,specfyc1,nfiles)
      call spectray(dc_dz,dc_dz,specep1,nfiles)
      call spectray(dcfluc,dcfluc,specep2,nfiles)
      call spectray(ducfluc,cfluc,spectp1,nfiles)
      call spectray(dwcfluc,cfluc,spectz1,nfiles)
      call spectrayi(cfluc,vcfluc,specsv1,nfiles)
      call spectray(fxfluc,dcfluc,specer1,nfiles)
      call spectray(fzfluc,dc_dz,specer2,nfiles)


      do i=0,imx
         do j=2,jmax-2,2
        specyuu(i,j/2)=specyuu(i,j/2)+(specuu1(i,j)+specuu1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyvv(i,j/2)=specyvv(i,j/2)+(specvv1(i,j)+specvv1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyww(i,j/2)=specyww(i,j/2)+(specww1(i,j)+specww1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyuw(i,j/2)=specyuw(i,j/2)+(specuw1(i,j)+specuw1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specypp(i,j/2)=specypp(i,j/2)+(specpp1(i,j)+specpp1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyqq(i,j/2)=specyqq(i,j/2)+(specqq1(i,j)+specqq1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specycc(i,j/2)=specycc(i,j/2)+(speccc1(i,j)+speccc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyec(i,j/2)=specyec(i,j/2)+(specec1(i,j)+specec1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyee(i,j/2)=specyee(i,j/2)+(specee1(i,j)+specee1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyGG(i,j/2)=specyGG(i,j/2)+(specGG1(i,j)+specGG1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyPe(i,j/2)=specyPe(i,j/2)+(specPe1(i,j)+specPe1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyPa(i,j/2)=specyPa(i,j/2)+(specPa1(i,j)+specPa1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyuc(i,j/2)=specyuc(i,j/2)+(specuc1(i,j)+specuc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specywc(i,j/2)=specywc(i,j/2)+(specwc1(i,j)+specwc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyqc(i,j/2)=specyqc(i,j/2)+(specqc1(i,j)+specqc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specykk(i,j/2)=specykk(i,j/2)+(speckk1(i,j)+speckk1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyfxc(i,j/2)=specyfxc(i,j/2)+(specfxc1(i,j)+specfxc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyfyc(i,j/2)=specyfyc(i,j/2)+(specfyc1(i,j)+specfyc1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specysv(i,j/2)=specysv(i,j/2)+(specsv1(i,j)+specsv1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specytp(i,j/2)=specytp(i,j/2)+(spectp1(i,j)+spectp1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specytz(i,j/2)=specytz(i,j/2)+(spectz1(i,j)+spectz1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyep1(i,j/2)=specyep1(i,j/2)+(specep1(i,j)+specep1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyep2(i,j/2)=specyep2(i,j/2)+(specep2(i,j)+specep2(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyer1(i,j/2)=specyer1(i,j/2)+(specer1(i,j)+specer1(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyer2(i,j/2)=specyer2(i,j/2)+(specer2(i,j)+specer2(i,j+1))/(jmax*nfiles)*2./(kmax)

       enddo

        specyuu(i,0)=specyuu(i,0)+(specuu1(i,1))/(jmax*kmax*nfiles)
        specyvv(i,0)=specyvv(i,0)+(specvv1(i,1))/(jmax*kmax*nfiles)
        specyww(i,0)=specyww(i,0)+(specww1(i,1))/(jmax*kmax*nfiles)
        specyuw(i,0)=specyuw(i,0)+(specuw1(i,1))/(jmax*kmax*nfiles)
        specypp(i,0)=specypp(i,0)+(specpp1(i,1))/(jmax*kmax*nfiles)
        specyqq(i,0)=specyqq(i,0)+(specqq1(i,1))/(jmax*nfiles*kmax)
        specycc(i,0)=specycc(i,0)+(speccc1(i,1))/(jmax*nfiles*kmax)
        specyec(i,0)=specyec(i,0)+(specec1(i,1))/(jmax*nfiles*kmax)
        specyee(i,0)=specyee(i,0)+(specee1(i,1))/(jmax*nfiles*kmax)
        specyGG(i,0)=specyGG(i,0)+(specGG1(i,1))/(jmax*nfiles*kmax)
        specyPe(i,0)=specyPe(i,0)+(specPe1(i,1))/(jmax*nfiles*kmax)
        specyPa(i,0)=specyPa(i,0)+(specPa1(i,1))/(jmax*nfiles*kmax)
        specyuc(i,0)=specyuc(i,0)+(specuc1(i,1))/(jmax*nfiles*kmax)
        specywc(i,0)=specywc(i,0)+(specwc1(i,1))/(jmax*nfiles*kmax)
        specyqc(i,0)=specyqc(i,0)+(specqc1(i,1))/(jmax*nfiles*kmax)
        specykk(i,0)=specykk(i,0)+(speckk1(i,1))/(jmax*nfiles*kmax)
        specyfxc(i,0)=specyfxc(i,0)+(specfxc1(i,1))/(jmax*nfiles*kmax)
        specyfyc(i,0)=specyfyc(i,0)+(specfyc1(i,1))/(jmax*nfiles)/(kmax)
        specysv(i,0)=specysv(i,0)+(specsv1(i,1))/(jmax*nfiles)/(kmax)
        specytp(i,0)=specytp(i,0)+(spectp1(i,1))/(jmax*nfiles)/(kmax)
        specytz(i,0)=specytz(i,0)+(spectz1(i,1))/(jmax*nfiles)/(kmax)
        specyep1(i,0)=specyep1(i,0)+(specep1(i,1))/(jmax*nfiles)/(kmax)
        specyep2(i,0)=specyep2(i,0)+(specep2(i,1))/(jmax*nfiles)/(kmax)
        specyer1(i,0)=specyer1(i,0)+(specer1(i,1))/(jmax*nfiles)/(kmax)
        specyer2(i,0)=specyer2(i,0)+(specer2(i,1))/(jmax*nfiles)/(kmax)

        specyuu(i,jmax/2)=specyuu(i,jmax/2)+(specuu1(i,jmax))/(jmax*kmax*nfiles)*2.
        specyvv(i,jmax/2)=specyvv(i,jmax/2)+(specvv1(i,jmax))/(jmax*kmax*nfiles)*2.
        specyww(i,jmax/2)=specyww(i,jmax/2)+(specww1(i,jmax))/(jmax*kmax*nfiles)*2.
        specyuw(i,jmax/2)=specyuw(i,jmax/2)+(specuw1(i,jmax))/(jmax*kmax*nfiles)*2.
        specypp(i,jmax/2)=specypp(i,jmax/2)+(specpp1(i,jmax))/(jmax*kmax*nfiles)*2.
        specyqq(i,jmax/2)=specyqq(i,jmax/2)+(specqq1(i,jmax))/(jmax*nfiles*kmax)*2.
        specycc(i,jmax/2)=specycc(i,jmax/2)+(speccc1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyec(i,jmax/2)=specyec(i,jmax/2)+(specec1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyee(i,jmax/2)=specyee(i,jmax/2)+(specee1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyGG(i,jmax/2)=specyGG(i,jmax/2)+(specGG1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyPe(i,jmax/2)=specyPe(i,jmax/2)+(specPe1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyPa(i,jmax/2)=specyPa(i,jmax/2)+(specPa1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyuc(i,jmax/2)=specyuc(i,jmax/2)+(specuc1(i,jmax))/(jmax*nfiles*kmax)*2.
        specywc(i,jmax/2)=specywc(i,jmax/2)+(specwc1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyqc(i,jmax/2)=specyqc(i,jmax/2)+(specqc1(i,jmax))/(jmax*nfiles*kmax)*2.
        specykk(i,jmax/2)=specykk(i,jmax/2)+(speckk1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyfxc(i,jmax/2)=specyfxc(i,jmax/2)+(specfxc1(i,jmax))/(jmax*nfiles*kmax)*2.
        specyfyc(i,jmax/2)=specyfyc(i,jmax/2)+(specfyc1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specysv(i,jmax/2)=specysv(i,jmax/2)+(specsv1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specytp(i,jmax/2)=specytp(i,jmax/2)+(spectp1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specytz(i,jmax/2)=specytz(i,jmax/2)+(spectz1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specyep1(i,jmax/2)=specyep1(i,jmax/2)+(specep1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specyep2(i,jmax/2)=specyep2(i,jmax/2)+(specep2(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specyer1(i,jmax/2)=specyer1(i,jmax/2)+(specer1(i,jmax))/(jmax*nfiles)/(kmax)*2.
        specyer2(i,jmax/2)=specyer2(i,jmax/2)+(specer2(i,jmax))/(jmax*nfiles)/(kmax)*2.

       enddo

      do i=1,imax
       yint_s(i,:)=specycc(i,:)
      enddo
      call deriv_s_c_m(imax,jmax/2+1,yint_s,yint_c,dr)
      do j=0,jmax/2
       yint_c(:,j)=yint_c(:,j)*mr_c(:)
      enddo
      call deriv_c_s_m(imax,jmax/2+1,yint_c,yder_s,dr)
      do j=0,jmax/2
       yder_s(:,j)=yder_s(:,j)*mr_s(:)
      enddo
      yint_s=0.
      yint_c=0.
      yder2_s=0.
      do i=1,imax
       yint_s(i,:)=specyfxc(i,:)
      enddo
      call deriv_s_c_m(imax,jmax/2+1,yint_s,yint_c,dr)
       do j=0,jmax/2
        yint_c(:,j)=yint_c(:,j)*mr_c(:)
       enddo
      call inter_c_s_m(imax,jmax/2+1,yint_c,yder2_s,dr)

      do i=1,imx
       specbudy(i,:,1)=-2*specyuc(i,:)*dcc(i)
       do j=0,jmax/2
        specbudy(i,j,2)=(specysv(i,j)*(j*2./1.5))
       enddo
       specbudy(i,:,3)=specytp(i,:)
       specbudy(i,:,4)=specytz(i,:)
       do j=0,jmax/2
        specbudy(i,j,5)=-specycc(i,j)*((j*2./1.5)**2.)*2./(Re*Pr)
       enddo
       specbudy(i,:,6)=-(specyep1(i,:))*2./(Re*Pr)
       specbudy(i,:,7)=-(specyep2(i,:))*2./(Re*Pr)
       specbudy(i,:,8)=(yder_s(i,:))/(Re*Pr)
       specbudy(i,:,9)=-(yder2_s(i,:))*2./(Re*Pr*Pl)
       specbudy(i,:,10)=(specyer1(i,:))*2./(Re*Pr*Pl)
       specbudy(i,:,11)=(specyer2(i,:))*2./(Re*Pr*Pl)
       do j=0,jmax/2
        specbudy(i,j,12)=(specyfyc(i,j)*(j*2./1.5))*2./(Re*Pr*Pl)
       enddo
       specbudy(i,:,13)=-(specyPe(i,:))*2./(Re*Pr*Pl)
       specbudy(i,:,14)=(specyPa(i,:))*2./(Re*Pr*Pl)
      enddo

      return
      end

      subroutine spectraz(p1x,p2x,spectra,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real spectra(0:imx,kmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wk(2*kmax+15),fz(kmax),fz2(kmax),dfz(kmax),dfz2(kmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,kmax/p_col)

      integer ier,nfiles

      call vrffti(kmax,wk)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)

      call transpose_x_to_y(p2x,p2y)
      call transpose_y_to_z(p2y,p2z)


       do i=0,imx
        do j=1,jmax/p_col
         do k=1,kmax
          fz(k) =p1z(i,j,k)
          fz2(k)=p2z(i,j,k)
         enddo
         call vrfftf(1,kmax,fz ,dfz ,1,wk)
         call vrfftf(1,kmax,fz2,dfz2,1,wk)
         do k=1,kmax
          p1z(i,j,k)=fz(k)*fz2(k)
        enddo
       enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         tmp(i,k)=tmp(i,k)+p1y(i,j,k)
        enddo
       enddo
      enddo
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         p1y(i,j,k)=tmp(i,k)
        enddo
       enddo
      enddo 

      call transpose_y_to_z(p1y,p1z)
      do k=1,kmax
       spectra(:,k) = p1z(:,1,k)
      end do


      return
      end



      subroutine spectray(p1x,p2x,spectra,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real spectra(0:imx,jmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wt(2*jmax+15),ft(jmax),ft2(jmax),dft(jmax),dft2(jmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,jmax/p_col)

      integer ier,nfiles

      call vrffti(jmax,wt)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_x_to_y(p2x,p2y)


       do i=0,imx
        do k=1,kmax/p_col
         do j=1,jmax
          ft(j) =p1y(i,j,k)
          ft2(j)=p2y(i,j,k)
         enddo
         call vrfftf(1,jmax,ft ,dft ,1,wt)
         call vrfftf(1,jmax,ft2,dft2,1,wt)
         do j=1,jmax
          p1y(i,j,k)=ft(j)*ft2(j)
        enddo
       enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         tmp(i,j)=tmp(i,j)+p1z(i,j,k)
        enddo
       enddo
      enddo
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         p1z(i,j,k)=tmp(i,j)
        enddo
       enddo
      enddo 

      call transpose_z_to_y(p1z,p1y)
      do j=1,jmax
       spectra(:,j) = p1y(:,j,1)
      end do


      return
      end

      subroutine spectrazi(p1x,p2x,spectra,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real spectra(0:imx,kmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wk(2*kmax+15),fz(kmax),fz2(kmax),dfz(kmax),dfz2(kmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,kmax/p_col)

      integer ier,nfiles

      call vrffti(kmax,wk)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)

      call transpose_x_to_y(p2x,p2y)
      call transpose_y_to_z(p2y,p2z)


       do i=0,imx
        do j=1,jmax/p_col
         do k=1,kmax
          fz(k) =p1z(i,j,k)
          fz2(k)=p2z(i,j,k)
         enddo
         call vrfftf(1,kmax,fz ,dfz ,1,wk)
         call vrfftf(1,kmax,fz2,dfz2,1,wk)
         p1z(i,j,1)=0
         p1z(i,j,kmax)=0
         do k=2,kmax-1
          if(mod(k,2).eq.0) then
           p1z(i,j,k)=-fz(k)*fz2(k+1)
          else
           p1z(i,j,k)=fz(k)*fz2(k-1)
          endif
        enddo
       enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         tmp(i,k)=tmp(i,k)+p1y(i,j,k)
        enddo
       enddo
      enddo
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         p1y(i,j,k)=tmp(i,k)
        enddo
       enddo
      enddo 

      call transpose_y_to_z(p1y,p1z)
      do k=1,kmax
       spectra(:,k) = p1z(:,1,k)
      end do


      return
      end



      subroutine spectrayi(p1x,p2x,spectra,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real spectra(0:imx,jmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wt(2*jmax+15),ft(jmax),ft2(jmax),dft(jmax),dft2(jmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,jmax/p_col)

      integer ier,nfiles

      call vrffti(jmax,wt)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_x_to_y(p2x,p2y)


       do i=0,imx
        do k=1,kmax/p_col
         do j=1,jmax
          ft(j) =p1y(i,j,k)
          ft2(j)=p2y(i,j,k)
         enddo
         call vrfftf(1,jmax,ft ,dft ,1,wt)
         call vrfftf(1,jmax,ft2,dft2,1,wt)
           p1y(i,1,k)=0.
           p1y(i,jmax,k)=0.
         do j=2,jmax-1
          if(mod(j,2).eq.0) then
           p1y(i,j,k)=-ft(j)*ft2(j+1)
          else
           p1y(i,j,k)=ft(j)*ft2(j-1)
          endif
        enddo
       enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         tmp(i,j)=tmp(i,j)+p1z(i,j,k)
        enddo
       enddo
      enddo
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         p1z(i,j,k)=tmp(i,j)
        enddo
       enddo
      enddo 

      call transpose_z_to_y(p1z,p1y)
      do j=1,jmax
       spectra(:,j) = p1y(:,j,1)
      end do


      return
      end

      subroutine yzderiv(pix,dpix_dt,dpix_dz)
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'

      real pix(0:imax+1,jmax/p_row,kmax/p_col)
      real dpix_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dpix_dz(0:imax+1,jmax/p_row,kmax/p_col)

      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax)
      real tmp (0:imx,jmax,kmax/p_col),tmpf (0:imx,jmax,kmax/p_col)
      real tmp2(0:imx,jmax,kmax/p_col),tmp2f(0:imx,jmax,kmax/p_col)
      real tmp3(0:imx,jmax/p_col,kmax),tmp3f(0:imx,jmax/p_col,kmax)
      integer ier
 
! -------------------------------------------------------------------

      call transpose_x_to_y(pix,tmp)

       do i=0,imx
         do k=1,kmax/p_col
          do j=1,jmax
            ft (j)=tmp (i,j,k)
          enddo
          call four1(jmax,ft,dft,dtheta)
          do j=1,jmax
           tmp (i,j,k)=dft (j)
          enddo
         enddo
      enddo
      call transpose_y_to_x(tmp,dpix_dt)
!-------------------------------------------------------------------

      call transpose_x_to_y(pix,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
      do i=0,imx
         do j=1,jmax/p_col
           do k=1,kmax
             fz (k)=tmp3 (i,j,k)
           enddo
           call four1(kmax,fz ,dfz ,dz)
           do k=1,kmax
            tmp3 (i,j,k)=dfz(k)
           enddo
          enddo
        enddo
        call transpose_z_to_y(tmp3,tmp2)
        call transpose_y_to_x(tmp2,dpix_dz)
!------------------------------------------------------------------
!-------------------------------------------------------------------

        return
        end

      subroutine corry(p1x,p2x,Ry,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real Ry(0:imx,jmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wt(2*jmax+15),ft(jmax),ft2(jmax),dft(jmax),dft2(jmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,jmax/p_col)

      integer ier,nfiles

      call vrffti(jmax,wt)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_x_to_y(p2x,p2y)


       do i=0,imx
        do k=1,kmax/p_col
         do j=1,jmax
          ft(j) =p1y(i,j,k)
          ft2(j)=p2y(i,j,k)
         enddo
         call vrfftf(1,jmax,ft ,dft ,1,wt)
         call vrfftf(1,jmax,ft2,dft2,1,wt)
           p1y(i,1,k)=ft(1)*ft2(1)
           p1y(i,jmax,k)=ft(jmax)*ft2(jmax)
         do j=2,jmax-1
          if(mod(j,2).eq.0) then
           p1y(i,j,k)=ft(j)*ft2(j)+ft(j+1)*ft2(j+1)
          else
           p1y(i,j,k)=ft(j)*ft2(j-1)-ft(j-1)*ft2(j)
          endif
        enddo
       enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         tmp(i,j)=tmp(i,j)+p1z(i,j,k)
        enddo
       enddo
      enddo
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         p1z(i,j,k)=tmp(i,j)
        enddo
       enddo
      enddo 

      call transpose_z_to_y(p1z,p1y)
      do i=0,imx
       do j=1,jmax
        ft(j)=p1y(i,j,1)
       enddo
       call vrfftb(1,jmax,ft,dft,1,wt)
       do j=1,jmax
        Ry(i,j)=ft(j)
       enddo       
      enddo


      return
      end

      subroutine corrz(p1x,p2x,Rz,nfiles)
      
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real Rz(0:imx,kmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,jmax,kmax/p_col)
      real p1z(0:imx,jmax/p_col,kmax)
      real wk(2*kmax+15),fz(kmax),fz2(kmax),dfz(kmax),dfz2(kmax)
      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
      real p2y(0:imx,jmax,kmax/p_col)
      real p2z(0:imx,jmax/p_col,kmax)
      real tmp(0:imx,kmax/p_col)

      integer ier,nfiles

      call vrffti(kmax,wk)

      tmp = 0 

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)

      call transpose_x_to_y(p2x,p2y)
      call transpose_y_to_z(p2y,p2z)


       do i=0,imx
        do j=1,jmax/p_col
         do k=1,kmax
          fz(k) =p1z(i,j,k)
          fz2(k)=p2z(i,j,k)
         enddo
         call vrfftf(1,kmax,fz ,dfz ,1,wk)
         call vrfftf(1,kmax,fz2,dfz2,1,wk)
         p1z(i,j,1)=fz(1)*fz2(1)
         p1z(i,j,kmax)=fz(kmax)*fz2(kmax)
         do k=2,kmax-1
          if(mod(k,2).eq.0) then
           p1z(i,j,k)=fz(k)*fz2(k)+fz(k+1)*fz2(k+1)
          else
           p1z(i,j,k)=fz(k)*fz2(k-1)-fz(k-1)*fz2(k)
          endif
        enddo
       enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         tmp(i,k)=tmp(i,k)+p1y(i,j,k)
        enddo
       enddo
      enddo
      do k=1,kmax/p_col
       do i=0,imx
        do j=1,jmax
         p1y(i,j,k)=tmp(i,k)
        enddo
       enddo
      enddo 

      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
       do k=1,kmax
        fz(k)=p1z(i,1,k)
       enddo
       call vrfftb(1,kmax,fz,dfz,1,wk)
       do k=1,kmax
        Rz(i,k) = fz(k)
       enddo
      enddo


      return
      end

      subroutine stateIG()
      implicit none
      include 'par_post.txt'
      include 'common.txt'

       do k=1,kmax/p_col
          do j=1,jmax/p_row
            do i=0,imax+1
              rhonew(i,j,k)=Tplus/(Tplus+cnew(i,j,k))
            enddo
          enddo
        enddo
      return
      end

      subroutine meanenergy(nfiles,rm,umf,cm,cmf,uflucf,cflucf,fxm,energym)
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      real uflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real cflucf(0:imax+1,jmax/p_row, kmax/p_col)
      real umf(0:imax),cmf(imax),cm(imax),fxm(imax),rm(imax)
      real energym(imax,4)
      real ufluc_c (0:imax,jmax/p_row, kmax/p_col)
      real ufluc_s (imax,jmax/p_row, kmax/p_col)
      real umean_c (0:imax)
      real umean_s (imax)
      real dc_c(0:imax),dc_s(imax)
      integer nfiles

      ! write(*,*) 'Inside job: -> ',nfiles
      ! terms in the mean integrated energy equation, as follows:
      !
      ! energym(i,1),  advection:     -->  rm * umf * cmf
      ! energym(i,2),  conduction:    -->  - d cm / dx * 1/RePr
      ! energym(i,3),  turb heat:     -->  rm * tilde(urf * crf)
      ! energym(i,4),  radiation:     -->  fxm* 1/RePrPl
 
      umean_c = umf
      do i=0,imax
       ufluc_c(i,:,:) = uflucf(i,:,:)
      enddo
      dc_s = cm

      call inter_c_s_m(imax,1,umean_c,umean_s,dr)
      call inter_c_s_m(imax,jmax*kmax/px,ufluc_c,ufluc_s,dr)
      call der1c_s_c6_m(imax,1,dc_s,dc_c,dr,Twb,Twt)
      dc_c = dc_c*mr_c
      call inter_c_s_m(imax,1,dc_c,dc_s,dr)

      do i=1,imax
       energym(i,1) = energym(i,1) + rm(i) * umf(i) * cmf(i) / (nfiles) 
       energym(i,2) = energym(i,2) - dc_s(i) / (nfiles*Re*Pr)
       energym(i,4) = energym(i,4) + fxm(i)  / (nfiles*Re*Pr*Pl)
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         energym(i,3) = energym(i,3) + rm(i) * ufluc_s(i,j,k)*cflucf(i,j,k) / (nfiles*jmax*kmax)
        enddo
       enddo
      enddo

      end


      subroutine tempvarbal(nfiles,rm,umf,cm,cmf,utmpf,ctmp,ctmpf,qm,qtmp,tempvar)
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      real utmpf(0:imax+1,jmax/p_row, kmax/p_col)
      real ctmpf(0:imax+1,jmax/p_row, kmax/p_col)
      real ctmp(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc_dt(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc_dz(0:imax+1,jmax/p_row, kmax/p_col)
      real dcfluc_c(0:imax,jmax/p_row, kmax/p_col)
      real qtmp(0:imax+1,jmax/p_row, kmax/p_col)
      real umf(0:imax),cmf(imax),cm(imax),rm(imax),qm(imax)
      real tempvar(imax,8), tempdum(imax,jmax/p_row,kmax/p_col,8)
      real ufluc_c (0:imax,jmax/p_row, kmax/p_col)
      real ufluc_s (imax,jmax/p_row, kmax/p_col)
      real umean_c (0:imax)
      real umean_s (imax)
      real dc_c(0:imax),dc_s(imax),ddc_s(imax)
      real dcf_c(0:imax),dcf_s(imax),hi1_s(imax),hi1_c(0:imax)
      real hiu_s(imax,jmax/p_row, kmax/p_col)
      real hiu_c(0:imax,jmax/p_row, kmax/p_col)
      real hic_s(imax,jmax/p_row, kmax/p_col)
      real hic_c(0:imax,jmax/p_row, kmax/p_col)
      real hi_s(imax,jmax/p_row, kmax/p_col)
      real hi_c(0:imax,jmax/p_row, kmax/p_col)
      real avg_crf(imax)      
      integer nfiles 

      ! terms in the temperature variance budget equation, as follows:
      !
      ! tempvar(i,1),  advection:     -->  - 1/2 rm * umf * d tilde(crf^2) / dx
      ! tempvar(i,2),  turb transp:   -->  - 1/2 d (rm * tilde(urf * crf^2) / dx
      ! tempvar(i,3),  production:    -->  - rm * tilde(urf * crf) d cmf / dx
      ! tempvar(i,4),  molec 1 :      -->  - d^2 cm / d x^2 * avg(crf) / RePr
      ! tempvar(i,5),  mol diff:      -->  d avg((d cr / dx) * crf) / dx / RePr 
      ! tempvar(i,6),  mol diss:      -->  - avg( cr * d^2 cr / dx_j^2 ) / RePr
      ! tempvar(i,7),  radiation1:    -->  qm * avg(crf) / RePrPl
      ! tempvar(i,8),  radiation2:    -->  - avg( qr * crf ) /RePrPl
 
      tempdum=0

      umean_c = umf
      do i=0,imax
       hiu_c(i,:,:) = utmpf(i,:,:)
      enddo
      dc_s = cm
      avg_crf = cm - cmf

      ! CALCULATING NEEDED DERIVATIVES
      hi1_s=0
      hi1_c=0
      call inter_c_s_m(imax,1,umean_c,umean_s,dr)
      call inter_c_s_m(imax,jmax*kmax/px,hiu_c,hiu_s,dr)
      call inter_s_c_m(imax,1,dc_s,dc_c,dr)
      call deriv_c_s_m(imax,1,dc_c,dc_s,dr)
      dc_s = dc_s*mr_s
      hi1_s = dc_s
      call inter_s_c_m(imax,1,hi1_s,hi1_c,dr)
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      ddc_s = hi1_s*mr_s
      hi1_s = cmf
      call inter_s_c_m(imax,1,hi1_s,hi1_c,dr)
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      dcf_s = hi1_s*mr_s
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         hic_s(i,j,k) = ctmp(i,j,k)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hic_s,hic_c,dr)
      call deriv_c_s_m(imax,jmax*kmax/px,hic_c,hic_s,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         dcfluc(i,j,k)=hic_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo

      ! CALCULATING BUDGET TERMS
      do i=1,imax
       tempdum(i,:,:,4) = - ddc_s(i) * avg_crf(i) / (nfiles*jmax*kmax*Re*Pr) 
       tempdum(i,:,:,7) = + qm(i) * avg_crf(i) / (nfiles*jmax*kmax*Re*Pr*Pl) 
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         tempdum(i,j,k,1) = -0.5*rhonew(i,j,k)*ctmpf(i,j,k)**2/
     +                            rm(i)/(nfiles*jmax*kmax)
         tempdum(i,j,k,2) = -0.5*rhonew(i,j,k)*hiu_s(i,j,k)*
     +                       ctmpf(i,j,k)**2/(nfiles*jmax*kmax)
         tempdum(i,j,k,3) = -dcf_s(i)*rhonew(i,j,k)*hiu_s(i,j,k)*
     +                       ctmpf(i,j,k)/(nfiles*jmax*kmax)
         tempdum(i,j,k,5) = dcfluc(i,j,k)*ctmpf(i,j,k)/(nfiles*jmax*kmax*Re*Pr)
         tempdum(i,j,k,8) = -qtmp(i,j,k)*ctmpf(i,j,k)/(nfiles*jmax*kmax*Re*Pr*Pl)
        enddo
       enddo
      enddo
 
      ! DERIVING TERM 1 (ADVECTION), 2 (TURB TRANSP) AND 5 (MOLEC DIFFUSION)
      hi_s = 0
      hi_c = 0
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k) = tempdum(i,j,k,1)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         tempdum(i,j,k,1)=hi_s(i,j,k)*mr_s(i)*rm(i)*umf(i) 
        enddo
       enddo
      enddo
      hi_s = 0
      hi_c = 0
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k) = tempdum(i,j,k,2)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         tempdum(i,j,k,2)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo
      hi_s = 0
      hi_c = 0
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k) = tempdum(i,j,k,5)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         tempdum(i,j,k,5)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo
      
      ! CALCULATING DISSIPATION TERM TEMPDUM(i,6)
      call yzderiv(ctmp,dcfluc_dt,dcfluc_dz)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         tempdum(i,j,k,6) = - dcfluc(i,j,k)   *  dcfluc(i,j,k)
     +                      - dcfluc_dt(i,j,k)*  dcfluc_dt(i,j,k)
     +                      - dcfluc_dz(i,j,k)*  dcfluc_dz(i,j,k)
         tempdum(i,j,k,6) = tempdum(i,j,k,6) / (nfiles*jmax*kmax*Re*Pr) 
        enddo
       enddo
      enddo

      ! REYNOLDS AVERAGING THE RESULTS
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         do h=1,8
          tempvar(i,h) = tempvar(i,h) + 2*tempdum(i,j,k,h)
         enddo
        enddo
       enddo
      enddo
      
      end

      subroutine turbflubal(nfiles,rm,um,umf,cm,cmf,qm,pm,utmp,utmpf,ctmp,ctmpf,qtmp,
     +                ptmp,vtmpf,wtmpf,vtmp,wtmp,turbflu)
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real utmp(0:imax+1,jmax/p_row,kmax/p_col)
      real utmpf(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmpf(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmpf(0:imax+1,jmax/p_row,kmax/p_col)
      real ctmp(0:imax+1,jmax/p_row, kmax/p_col)
      real ctmpf(0:imax+1,jmax/p_row,kmax/p_col)
      real ptmp(0:imax+1,jmax/p_row,kmax/p_col)
      real qtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real rm(imax),um(0:imax),umf(0:imax),cm(imax),cmf(imax),qm(imax),pm(imax)
      real turbflu(imax,15)
      real turbdum(0:imax+1,jmax/p_row, kmax/p_col,15)
      real hi_s(imax    ,jmax/p_row, kmax/p_col)
      real hi_c(0:imax  ,jmax/p_row, kmax/p_col)
      real utmpf_s(0:imax+1,jmax/p_row,kmax/p_col)
      real utmp_s(0:imax+1,jmax/p_row,kmax/p_col)
      real dutmpf(0:imax+1,jmax/p_row,kmax/p_col)
      real dutmp(0:imax+1,jmax/p_row,kmax/p_col)
      real dvtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real dwtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real dctmp(0:imax+1,jmax/p_row,kmax/p_col)
      real du_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real du_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real hiu_s(imax    ,jmax/p_row,kmax/p_col)
      real hiv_s(imax    ,jmax/p_row,kmax/p_col)
      real hiw_s(imax    ,jmax/p_row,kmax/p_col)
      real hic_s(imax    ,jmax/p_row,kmax/p_col)
      real hiu_c(0:imax  ,jmax/p_row,kmax/p_col)
      real hiv_c(0:imax  ,jmax/p_row,kmax/p_col)
      real hiw_c(0:imax  ,jmax/p_row,kmax/p_col)
      real hic_c(0:imax  ,jmax/p_row,kmax/p_col)
      real hi1_s(imax),hi1_c(0:imax)
      real dumf(imax),dum(imax),dcmf(imax),trace
      real t11m(imax),avg_cf(imax),avg_uf(imax)
      real tens(0:imax+1,jmax/p_row,kmax/p_col,3,3)
      real tmp(0:imx,jmax,kmax/p_col)
      integer nfiles, ier

      ! terms in the turbulent heat flux budget equation, as follows:
      !
      ! turbflu(i,1),  advection:     -->  - rm * umf * d tilde(crf*urf) / dx
      ! turbflu(i,2),  production1:   -->  - rm * tilde(urf * crf) d umf / dx
      ! turbflu(i,3),  production2:   -->  - rm * tilde(urf^2) d cmf / dx
      ! turbflu(i,4),  turb transp:   -->  - d rmf * tilde(vrf^2 * crf) / dx
      ! turbflu(i,5),  visc diff:     -->  d avg( t11r * cr) / d x / Re
      ! turbflu(i,6),  visc diss:     -->  - avg( tijr * d cr / dx_j ) / Re
      ! turbflu(i,7),  viscous1:      -->  - d t11m / dx * avg(crf) / Re
      ! turbflu(i,8),  press diff:    -->  - d avg( pr * cr) / d x 
      ! turbflu(i,9),  press diss:    -->  avg( pr * d cr / d x )
      ! turbflu(i,10), pressure1:     -->  d pm / dx * avg(crf)
      ! turbflu(i,11), mol diff:      -->  d avg( Fcr * ur) / d x / RePr
      ! turbflu(i,12), mol diss:      -->  - avg( Fcr * d ur / d x_j ) / RePr
      ! turbflu(i,13), molecular1:    -->  - d Fcm / dx * avg(urf) / RePr
      ! turbflu(i,14), rad term:      -->  - avg( qr * ur ) / RePrPl
      ! turbflu(i,15), radiation1:    -->  qm * avg(urf) / RePrPl

      
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_c(i,j,k)=utmp(i,j,k)
        enddo
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         utmp_s(i,j,k)=hi_s(i,j,k)
        enddo
       enddo
      enddo
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_c(i,j,k)=utmpf(i,j,k)
        enddo
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         utmpf_s(i,j,k)=hi_s(i,j,k)
        enddo
       enddo
      enddo
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         dutmpf(i,j,k)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo
      hi1_c = umf
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      dumf = hi1_s*mr_s
      hi1_c = um
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      dum = hi1_s*mr_s
      hi1_s = cmf
      call inter_s_c_m(imax,1,hi1_s,hi1_c,dr)
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      dcmf = hi1_s*mr_s
      avg_cf = cm - cmf
      hi1_c  = um - umf
      call inter_c_s_m(imax,1,hi1_c,hi1_s,dr)
      avg_uf = hi1_s

      ! CALCULATING THE LHS (WITH THE SIGNS IN THE RHS)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,1) = - rhonew(i,j,k)*ctmpf(i,j,k)*utmpf_s(i,j,k) / rm(i)
     +                      /(nfiles*jmax*kmax)
         turbdum(i,j,k,2) = - rhonew(i,j,k)*ctmpf(i,j,k)*utmpf_s(i,j,k)*dumf(i)
     +                      /(nfiles*jmax*kmax)
         turbdum(i,j,k,3) = - rhonew(i,j,k)*utmpf_s(i,j,k)**2*dcmf(i)
     +                      /(nfiles*jmax*kmax)
         turbdum(i,j,k,4) = - rhonew(i,j,k)*utmpf_s(i,j,k)**2*ctmpf(i,j,k)
     +                      /(nfiles*jmax*kmax)
        enddo
       enddo
      enddo
      
      ! DERIVING TERMS 1 (ADVECTION) 4 (TURBULENT TRANSPORT)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k)=turbdum(i,j,k,1)*rm(i)*umf(i)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,1)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k)=turbdum(i,j,k,4)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,4)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo

      ! CALCULATING THE FLUCTUATING STRESS TENSOR AND t22 MEAN

      call yzderiv(utmp  ,du_dt,du_dz)
      call yzderiv(vtmp  ,dv_dt,dv_dz)
      call yzderiv(wtmp  ,dw_dt,dw_dz)
      call yzderiv(ctmp  ,dc_dt,dc_dz)
      do j=1,jmax/p_row
       do k=1,kmax/p_col
        do i=1,imax
         hiv_s(i,j,k)=vtmp(i,j,k)
         hiw_s(i,j,k)=wtmp(i,j,k)
         hic_s(i,j,k)=ctmp(i,j,k)
        enddo
        do i=0,imax
         hiu_c(i,j,k)=utmp(i,j,k)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hiv_s,hiv_c,dr) 
      call inter_s_c_m(imax,jmax*kmax/px,hiw_s,hiw_c,dr) 
      call inter_s_c_m(imax,jmax*kmax/px,hic_s,hic_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hiu_c,hiu_s,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hiv_c,hiv_s,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hiw_c,hiw_s,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hic_c,hic_s,dr) 
      do j=1,jmax/p_row
       do k=1,kmax/p_col
        do i=1,imax
         dutmp(i,j,k)=hiu_s(i,j,k)*mr_s(i)
         dvtmp(i,j,k)=hiv_s(i,j,k)*mr_s(i)
         dwtmp(i,j,k)=hiw_s(i,j,k)*mr_s(i)
         dctmp(i,j,k)=hic_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo
      do i=1,imax
       t11m(i) = 4./3. * dum(i) 
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         trace = 2./3. * (dutmp(i,j,k) + dv_dt(i,j,k) + dw_dz(i,j,k))
         tens(i,j,k,1,1)= 2*dutmp(i,j,k) - trace 
         tens(i,j,k,1,2)=   du_dt(i,j,k) + dvtmp(i,j,k)  
         tens(i,j,k,1,3)=   du_dz(i,j,k) + dwtmp(i,j,k) 
         tens(i,j,k,2,1)=   dvtmp(i,j,k) + du_dt(i,j,k) 
         tens(i,j,k,2,2)= 2*dv_dt(i,j,k) - trace 
         tens(i,j,k,2,3)=   dv_dz(i,j,k) + dw_dt(i,j,k) 
         tens(i,j,k,3,1)=   dwtmp(i,j,k) + du_dz(i,j,k) 
         tens(i,j,k,3,2)=   dw_dt(i,j,k) + dv_dz(i,j,k) 
         tens(i,j,k,3,3)= 2*dw_dz(i,j,k) - trace 
        enddo
       enddo
      enddo

      ! CALCULATION OF TERMS 5 TO 7 (VISCOUS TERMS)
      hi1_s = t11m 
      call inter_s_c_m(imax,1,hi1_s,hi1_c,dr)
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,5) = tens(i,j,k,1,1) * ctmp(i,j,k) 
     +                      /(nfiles*jmax*kmax*Re)
         turbdum(i,j,k,6) = ( - tens(i,j,k,1,1) * dctmp(i,j,k) 
     +                        - tens(i,j,k,1,2) * dc_dt(i,j,k) 
     +                        - tens(i,j,k,1,3) * dc_dz(i,j,k) ) 
     +                      /(nfiles*jmax*kmax*Re)
         turbdum(i,j,k,7) = - hi1_s(i)*mr_s(i) * avg_cf(i) 
     +                      /(nfiles*jmax*kmax*Re)
        enddo
       enddo
      enddo
      ! DERIVING TERM 5 (VISCOUS DIFFUSION)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k)=turbdum(i,j,k,5)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,5)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo

      ! CALCULATION OF TERMS 8 TO 10 (PRESSURE TERMS)
      hi1_s = pm 
      call inter_s_c_m(imax,1,hi1_s,hi1_c,dr)
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,8) = - ptmp(i,j,k) * ctmp(i,j,k) 
     +                      /(nfiles*jmax*kmax)
         turbdum(i,j,k,9) = ptmp(i,j,k) * dctmp(i,j,k)  
     +                      /(nfiles*jmax*kmax)
         turbdum(i,j,k,10) = - hi1_s(i)*mr_s(i) * avg_cf(i) 
     +                      /(nfiles*jmax*kmax)
        enddo
       enddo
      enddo
      ! DERIVING TERM 8 (PRESSURE DIFFUSION)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k)=turbdum(i,j,k,8)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,8)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo

      ! CALCULATION OF TERMS 11 TO 13 (MOLECULAR TERMS)
      hi1_s = cm 
      call deriv_s_c_m(imax,1,hi1_s,hi1_c,dr)
      hi1_c = hi1_c * mr_c
      call deriv_c_s_m(imax,1,hi1_c,hi1_s,dr)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,11) = - utmp_s(i,j,k) * dctmp(i,j,k) 
     +                      /(nfiles*jmax*kmax*Re*Pr)
         turbdum(i,j,k,12) = -(dctmp(i,j,k) * dutmp(i,j,k) +
     +                         dc_dt(i,j,k) * dv_dt(i,j,k) +
     +                         dc_dz(i,j,k) * dw_dz(i,j,k) )  
     +                      /(nfiles*jmax*kmax*Re*Pr)
         turbdum(i,j,k,13) = - hi1_s(i)*mr_s(i) * avg_uf(i) 
     +                      /(nfiles*jmax*kmax*Re*Pr)
        enddo
       enddo
      enddo
      ! DERIVING TERM 11 (MOLECULAR DIFFUSION)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         hi_s(i,j,k)=turbdum(i,j,k,11)
        enddo
       enddo
      enddo
      call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr) 
      call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,11)=hi_s(i,j,k)*mr_s(i)
        enddo
       enddo
      enddo

      ! CALCULATION OF TERMS 14 AND 15 (RADIATION TERMS) 
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         turbdum(i,j,k,14) = - qtmp(i,j,k) * utmp_s(i,j,k) 
     +                      /(nfiles*jmax*kmax*Re*Pr*Pl)
         turbdum(i,j,k,15) = qm(i) * avg_uf(i) 
     +                      /(nfiles*jmax*kmax*Re*Pr*Pl)
        enddo
       enddo
      enddo

      ! SUMMATION AND AVERAGING
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         do h=1,15
          turbflu(i,h) = turbflu(i,h) + turbdum(i,j,k,h)
         enddo 
        enddo
       enddo
      enddo

      end


      subroutine quad_analysis(Q,fluc1,fluc2,nfiles)

      implicit none
      include 'par_post.txt'
      include 'common.txt'
      real fluc1(0:imax+1,jmax/p_row, kmax/p_col)
      real fluc2(0:imax+1,jmax/p_row, kmax/p_col)
      real Q(imax,4)

      integer nfiles
         do k=1,kmax/p_col
          do j=1,jmax/p_row
            do i=1,imax

             if(fluc1(i,j,k).ge.0.and.fluc2(i,j,k).ge.0) Q(i,1)=Q(i,1)+(fluc1(i,j,k))*(fluc2(i,j,k))/((jmax)*(kmax)*nfiles)
             if(fluc1(i,j,k).ge.0.and.fluc2(i,j,k).lt.0) Q(i,2)=Q(i,2)+(fluc1(i,j,k))*(fluc2(i,j,k))/((jmax)*(kmax)*nfiles)
             if(fluc1(i,j,k).lt.0.and.fluc2(i,j,k).ge.0) Q(i,3)=Q(i,3)+(fluc1(i,j,k))*(fluc2(i,j,k))/((jmax)*(kmax)*nfiles)
             if(fluc1(i,j,k).lt.0.and.fluc2(i,j,k).lt.0) Q(i,4)=Q(i,4)+(fluc1(i,j,k))*(fluc2(i,j,k))/((jmax)*(kmax)*nfiles)

            enddo
          enddo
        enddo

      end
