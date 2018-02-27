     
      subroutine specbalyz(cfluc,ufluc_s,vfluc,wfluc,fxfluc,fyfluc,fzfluc,cm,specbudy,specbudz,
     +                    speczuc,speczcc,nfiles,rnk)

      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'
      real cfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dc_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dcfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real ufluc_s(0:imax+1,jmax/p_row,kmax/p_col)
      real ucfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real ducfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real fxfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real fyfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real fzfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real vfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real wfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real vcfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real wcfluc(0:imax+1,jmax/p_row,kmax/p_col)
      real cm(imax),dc(0:imax),dcc(imax)

      real specbudy(imax,0:jmax/2,9)
      real specbudz(imax,0:kmax/2,9)

      real specccy(0:imx,jmax)
      real specdcdz(0:imx,jmax)
      real specucy(0:imx,jmax)
      real specfxcy(0:imx,jmax)
      real specfycy(0:imx,jmax)
      real spectpy(0:imx,jmax)
      real specsvy(0:imx,jmax)
      real specepsy(0:imx,jmax)
      real specepsry(0:imx,jmax)
      real specepsr2y(0:imx,jmax)

      real specccz(0:imx,kmax)
      real specdcdy(0:imx,kmax)
      real specucz(0:imx,kmax)
      real specfxcz(0:imx,kmax)
      real specfzcz(0:imx,kmax)
      real spectpz(0:imx,kmax)
      real specswz(0:imx,kmax)
      real specepsz(0:imx,kmax)
      real specepsrz(0:imx,kmax)
      real specepsr2z(0:imx,kmax)

      real specycc(0:imx,0:jmax/2)
      real specydc(0:imx,0:jmax/2)
      real specyuc(0:imx,0:jmax/2)
      real specyfxc(0:imx,0:jmax/2)
      real specyfyc(0:imx,0:jmax/2)
      real specytp(0:imx,0:jmax/2)
      real specysv(0:imx,0:jmax/2)
      real specyeps(0:imx,0:jmax/2)
      real specyepsr(0:imx,0:jmax/2)
      real specyepsr2(0:imx,0:jmax/2)

      real speczcc(0:imx,0:kmax/2)
      real speczdc(0:imx,0:kmax/2)
      real speczuc(0:imx,0:kmax/2)
      real speczfxc(0:imx,0:kmax/2)
      real speczfzc(0:imx,0:kmax/2)
      real specztp(0:imx,0:kmax/2)
      real speczsw(0:imx,0:kmax/2)
      real speczeps(0:imx,0:kmax/2)
      real speczepsr(0:imx,0:kmax/2)
      real speczepsr2(0:imx,0:kmax/2)

      real interp_c(0:imax,jmax/p_row,kmax/p_col)
      real interp_s(imax,jmax/p_row,kmax/p_col)
      real deriv_s(imax,jmax/p_row,kmax/p_col)
      real yint_c(0:imax,0:jmax/2)
      real yint_s(imax,0:jmax/2)
      real yder_s(imax,0:jmax/2)
      real yder2_s(imax,0:jmax/2)
      real zint_c(0:imax,0:kmax/2)
      real zint_s(imax,0:kmax/2)
      real zder_s(imax,0:kmax/2)
      real zder2_s(imax,0:kmax/2)

      real stime
      integer nfiles,rnk,rs

      ducfluc=0
      dcfluc=0
      ucfluc=ufluc_s*cfluc
      vcfluc=wfluc*cfluc
      wcfluc=wfluc*cfluc
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

      call der1c_s_c6_m(imax,1,cm,dc,dr,Twb,Twt)
      dc=dc*mr_c
      call inter_c_s_m(imax,1,dc,dcc,dr)

      call yzderiv(cfluc,dc_dt,dc_dz)


      if(rnk.eq.0) write(*,*) 'entering 2D spectra'
      stime=MPI_WTIME()
      call spectray(cfluc,cfluc,specccy,nfiles) 
      call spectray(dc_dz,dc_dz,specdcdz,nfiles) 
      call spectray(ducfluc,cfluc,spectpy,nfiles) 
      call spectrayi(vcfluc,cfluc,specsvy,nfiles) 
      call spectray(ufluc_s,cfluc,specucy,nfiles) 
      call spectray(fxfluc,cfluc,specfxcy,nfiles) 
      call spectrayi(fyfluc,cfluc,specfycy,nfiles) 
      call spectray(dcfluc,dcfluc,specepsy,nfiles) 
      call spectray(fxfluc,dcfluc,specepsry,nfiles) 
      call spectray(fzfluc,dc_dz,specepsr2y,nfiles) 

      call spectraz(cfluc,cfluc,specccz,nfiles) 
      call spectraz(dc_dt,dc_dt,specdcdy,nfiles) 
      call spectraz(ducfluc,cfluc,spectpz,nfiles) 
      call spectrazi(wcfluc,cfluc,specswz,nfiles) 
      call spectraz(ufluc_s,cfluc,specucz,nfiles) 
      call spectraz(fxfluc,cfluc,specfxcz,nfiles) 
      call spectrazi(fzfluc,cfluc,specfzcz,nfiles) 
      call spectraz(dcfluc,dcfluc,specepsz,nfiles) 
      call spectraz(fxfluc,dcfluc,specepsrz,nfiles) 
      call spectraz(fyfluc,dc_dt,specepsr2z,nfiles) 

      if(rnk.eq.0) write(*,*) 'Time in spectra: ',-stime+MPI_WTIME()

      do i=0,imx
       do j=2,jmax-2,2
        specycc(i,j/2)=specycc(i,j/2)+(specccy(i,j)+specccy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specydc(i,j/2)=specydc(i,j/2)+(specdcdz(i,j)+specdcdz(i,j+1))/(jmax*nfiles)*2./(kmax)
        specytp(i,j/2)=specytp(i,j/2)+(spectpy(i,j)+spectpy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specysv(i,j/2)=specysv(i,j/2)+(specsvy(i,j)+specsvy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyuc(i,j/2)=specyuc(i,j/2)+(specucy(i,j)+specucy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyfxc(i,j/2)=specyfxc(i,j/2)+(specfxcy(i,j)+specfxcy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyfyc(i,j/2)=specyfyc(i,j/2)+(specfycy(i,j)+specfycy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyeps(i,j/2)=specyeps(i,j/2)+(specepsy(i,j)+specepsy(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyepsr(i,j/2)=specyepsr(i,j/2)+(specepsry(i,j)+specepsry(i,j+1))/(jmax*nfiles)*2./(kmax)
        specyepsr2(i,j/2)=specyepsr2(i,j/2)+(specepsr2y(i,j)+specepsr2y(i,j+1))/(jmax*nfiles)*2./(kmax)
       enddo

       specycc(i,0)=specycc(i,0)+(specccy(i,1))/(jmax*nfiles)*2./(kmax)
       specydc(i,0)=specydc(i,0)+(specdcdz(i,1))/(jmax*nfiles)*2./(kmax)
       specytp(i,0)=specytp(i,0)+(spectpy(i,1))/(jmax*nfiles)*2./(kmax)
       specysv(i,0)=specysv(i,0)+(specsvy(i,1))/(jmax*nfiles)*2./(kmax)
       specyuc(i,0)=specyuc(i,0)+(specucy(i,1))/(jmax*nfiles)*2./(kmax)
       specyfxc(i,0)=specyfxc(i,0)+(specfxcy(i,1))/(jmax*nfiles)*2./(kmax)
       specyfyc(i,0)=specyfyc(i,0)+(specfycy(i,1))/(jmax*nfiles)*2./(kmax)
       specyeps(i,0)=specyeps(i,0)+(specepsy(i,1))/(jmax*nfiles)*2./(kmax)
       specyepsr(i,0)=specyepsr(i,0)+(specepsry(i,1))/(jmax*nfiles)*2./(kmax)
       specyepsr2(i,0)=specyepsr2(i,0)+(specepsr2y(i,1))/(jmax*nfiles)*2./(kmax)
       
       specycc(i,jmax/2)=specycc(i,jmax/2)+(specccy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specydc(i,jmax/2)=specydc(i,jmax/2)+(specdcdz(i,jmax))/(jmax*nfiles)*2./(kmax)
       specytp(i,jmax/2)=specytp(i,jmax/2)+(spectpy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specysv(i,jmax/2)=specysv(i,jmax/2)+(specsvy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyuc(i,jmax/2)=specyuc(i,jmax/2)+(specucy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyfxc(i,jmax/2)=specyfxc(i,jmax/2)+(specfxcy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyfyc(i,jmax/2)=specyfyc(i,jmax/2)+(specfycy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyeps(i,jmax/2)=specyeps(i,jmax/2)+(specepsy(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyepsr(i,jmax/2)=specyepsr(i,jmax/2)+(specepsry(i,jmax))/(jmax*nfiles)*2./(kmax)
       specyepsr2(i,jmax/2)=specyepsr2(i,jmax/2)+(specepsr2y(i,jmax))/(jmax*nfiles)*2./(kmax)

      enddo

      do i=0,imx
       do k=2,kmax-2,2
        speczcc(i,k/2)=speczcc(i,k/2)+(specccz(i,k)+specccz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczdc(i,k/2)=speczdc(i,k/2)+(specdcdy(i,k)+specdcdy(i,k+1))/(jmax*nfiles)*2./(kmax)
        specztp(i,k/2)=specztp(i,k/2)+(spectpz(i,k)+spectpz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczsw(i,k/2)=speczsw(i,k/2)+(specswz(i,k)+specswz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczuc(i,k/2)=speczuc(i,k/2)+(specucz(i,k)+specucz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczfxc(i,k/2)=speczfxc(i,k/2)+(specfxcz(i,k)+specfxcz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczfzc(i,k/2)=speczfzc(i,k/2)+(specfzcz(i,k)+specfzcz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczeps(i,k/2)=speczeps(i,k/2)+(specepsz(i,k)+specepsz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczepsr(i,k/2)=speczepsr(i,k/2)+(specepsrz(i,k)+specepsrz(i,k+1))/(jmax*nfiles)*2./(kmax)
        speczepsr2(i,k/2)=speczepsr2(i,k/2)+(specepsr2z(i,k)+specepsr2z(i,k+1))/(jmax*nfiles)*2./(kmax)
       enddo

       speczcc(i,0)=speczcc(i,0)+(specccz(i,1))/(jmax*nfiles)*2./(kmax)
       speczdc(i,0)=speczdc(i,0)+(specdcdy(i,1))/(jmax*nfiles)*2./(kmax)
       specztp(i,0)=specztp(i,0)+(spectpz(i,1))/(jmax*nfiles)*2./(kmax)
       speczsw(i,0)=speczsw(i,0)+(specswz(i,1))/(jmax*nfiles)*2./(kmax)
       speczuc(i,0)=speczuc(i,0)+(specucz(i,1))/(jmax*nfiles)*2./(kmax)
       speczfxc(i,0)=speczfxc(i,0)+(specfxcz(i,1))/(jmax*nfiles)*2./(kmax)
       speczfzc(i,0)=speczfzc(i,0)+(specfzcz(i,1))/(jmax*nfiles)*2./(kmax)
       speczeps(i,0)=speczeps(i,0)+(specepsz(i,1))/(jmax*nfiles)*2./(kmax)
       speczepsr(i,0)=speczepsr(i,0)+(specepsrz(i,1))/(jmax*nfiles)*2./(kmax)
       speczepsr2(i,0)=speczepsr2(i,0)+(specepsr2z(i,1))/(jmax*nfiles)*2./(kmax)

       speczcc(i,kmax/2)=speczcc(i,kmax/2)+(specccz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczdc(i,kmax/2)=speczdc(i,kmax/2)+(specdcdy(i,kmax))/(jmax*nfiles)*2./(kmax)
       specztp(i,kmax/2)=specztp(i,kmax/2)+(spectpz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczsw(i,kmax/2)=speczsw(i,kmax/2)+(specswz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczuc(i,kmax/2)=speczuc(i,kmax/2)+(specucz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczfxc(i,kmax/2)=speczfxc(i,kmax/2)+(specfxcz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczfzc(i,kmax/2)=speczfzc(i,kmax/2)+(specfzcz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczeps(i,kmax/2)=speczeps(i,kmax/2)+(specepsz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczepsr(i,kmax/2)=speczepsr(i,kmax/2)+(specepsrz(i,kmax))/(jmax*nfiles)*2./(kmax)
       speczepsr2(i,kmax/2)=speczepsr2(i,kmax/2)+(specepsr2z(i,kmax))/(jmax*nfiles)*2./(kmax)

      enddo

      do i=1,imax
       yint_s(i,:)=specycc(i,:)
       zint_s(i,:)=speczcc(i,:)
      enddo
      call deriv_s_c_m(imax,jmax/2+1,yint_s,yint_c,dr)
      call deriv_s_c_m(imax,kmax/2+1,zint_s,zint_c,dr)
      do j=0,jmax/2
       yint_c(:,j)=yint_c(:,j)*mr_c(:)
      enddo
      do k=0,kmax/2
       zint_c(:,k)=zint_c(:,k)*mr_c(:)
      enddo
      call deriv_c_s_m(imax,jmax/2+1,yint_c,yder_s,dr)
      call deriv_c_s_m(imax,kmax/2+1,zint_c,zder_s,dr)
      do j=0,jmax/2
       yder_s(:,j)=yder_s(:,j)*mr_s(:)
      enddo
      do k=0,kmax/2
       zder_s(:,k)=zder_s(:,k)*mr_s(:)
      enddo
      do i=1,imax
       yint_s(i,:)=specyfxc(i,:)
       zint_s(i,:)=speczfxc(i,:)
      enddo
      call deriv_s_c_m(imax,jmax/2+1,yint_s,yint_c,dr)
      call deriv_s_c_m(imax,kmax/2+1,zint_s,zint_c,dr)
       do j=1,jmax/2
        yint_c(:,j)=yint_c(:,j)*mr_c(:)
       enddo
       do k=1,kmax/2
        zint_c(:,k)=zint_c(:,k)*mr_c(:)
       enddo
      call inter_c_s_m(imax,jmax/2+1,yint_c,yder2_s,dr)
      call inter_c_s_m(imax,kmax/2+1,zint_c,zder2_s,dr)

      if(rnk.eq.0) write(*,*) 'writing spectral budgets'
      do i=1,imax
       specbudy(i,:,1)=-2*specyuc(i,:)*dcc(i)
       do j=0,jmax/2
        specbudy(i,j,2)=+(specysv(i,j)*(-j*2./1.5))
       enddo
       specbudy(i,:,3)=-specytp(i,:)
       do j=0,jmax/2
        specbudy(i,j,4)=- specycc(i,j)*((j*2./1.5)**2.)/(Re*Pr) 
       enddo
       specbudy(i,:,5)=- (specyeps(i,:)+specydc(i,:))*2./(Re*Pr)
       specbudy(i,:,6)=+ (yder_s(i,:))/(Re*Pr)
       specbudy(i,:,7)=- (yder2_s(i,:))*2/(Re*Pr*Pl)
       specbudy(i,:,8)=+ (specyepsr(i,:)+specyepsr2(i,:))*2/(Re*Pr*Pl)
       do j=0,jmax/2
        specbudy(i,j,9)=+ (specyfyc(i,j)*(-j*2./1.5))*2/(Re*Pr*Pl)
       enddo
      enddo

      do i=1,imax
       specbudz(i,:,1)= - speczuc(i,:)*dcc(i)
       do k=0,kmax/2
        specbudz(i,k,2)=+ (speczsw(i,k)*(-k*2./4.))
       enddo
       specbudz(i,:,3)=- specztp(i,:)
       do k=0,kmax/2
        specbudz(i,k,4)=- speczcc(i,k)*((k*2./4.)**2.)/(Re*Pr) 
       enddo
       specbudz(i,:,5)=- (speczeps(i,:)+speczdc(i,:))*2./(Re*Pr)
       specbudz(i,:,6)=+ (zder_s(i,:))/(Re*Pr)
       specbudz(i,:,7)=- (zder2_s(i,:))*2/(Re*Pr*Pl)
       specbudz(i,:,8)=+ (speczepsr(i,:)+speczepsr2(i,:))*2/(Re*Pr*Pl)
       do k=0,kmax/2
        specbudz(i,k,9)=+ (speczfzc(i,k)*(-k*2./4.))*2/(Re*Pr*Pl)
       enddo
      enddo

      return
      end


 
!      subroutine specbal(cfluc,ufluc_s,vfluc,wfluc,fxfluc,fyfluc,fzfluc,qfluc,efluc,cm,specbud,specnorm,speccc,
!     +                   specccint,specqqint,speceeint,specqcint,specww,nfiles,rnk)
!
!      use decomp_2d
!      implicit none
!      include 'par_post.txt'
!      include 'common.txt'
!      include 'mpif.h'
!      real cfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real dcfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real ufluc_s(0:imax+1,jmax/p_row,kmax/p_col)
!      real ucfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real ducfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real qfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real fxfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real fyfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real fzfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real vfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real wfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real efluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real vcfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real wcfluc(0:imax+1,jmax/p_row,kmax/p_col)
!      real cm(imax),dc(0:imax),dcc(imax)
!
!      real specbud(imax,hmax,9)
!      real specnorm(imax,9)
!
!      real speccc(0:imx,jmax/2,kmax/2)
!      real specqq(0:imx,jmax/2,kmax/2)
!      real specee(0:imx,jmax/2,kmax/2)
!      real specqc(0:imx,jmax/2,kmax/2)
!      real speccci(0:imx,jmax/2,kmax/2)
!      real specww(0:imx,jmax/2,kmax/2)
!      real specwwi(0:imx,jmax/2,kmax/2)
!      real specuc(0:imx,jmax/2,kmax/2)
!      real specuci(0:imx,jmax/2,kmax/2)
!      real specfxc(0:imx,jmax/2,kmax/2)
!      real specfxci(0:imx,jmax/2,kmax/2)
!      real specfyc(0:imx,jmax/2,kmax/2)
!      real specfyci(0:imx,jmax/2,kmax/2)
!      real specfzc(0:imx,jmax/2,kmax/2)
!      real specfzci(0:imx,jmax/2,kmax/2)
!      real spectp(0:imx,jmax/2,kmax/2)
!      real spectpi(0:imx,jmax/2,kmax/2)
!      real spectsv(0:imx,jmax/2,kmax/2)
!      real spectsvi(0:imx,jmax/2,kmax/2)
!      real spectsw(0:imx,jmax/2,kmax/2)
!      real spectswi(0:imx,jmax/2,kmax/2)
!      real speceps(0:imx,jmax/2,kmax/2)
!      real specepsi(0:imx,jmax/2,kmax/2)
!      real specepsr(0:imx,jmax/2,kmax/2)
!      real specepsri(0:imx,jmax/2,kmax/2)
!
!      real interp_c(0:imax,jmax/p_row,kmax/p_col)
!      real interp_s(imax,jmax/p_row,kmax/p_col)
!      real deriv_s(imax,jmax/p_row,kmax/p_col)
!      real int_c(0:imax,hmax)
!      real int_s(imax,hmax)
!      real der_s(imax,hmax)
!      real der2_s(imax,hmax)
!      real dumint(0:imx,jmax/2)
!
!      real specccint(0:imx,hmax)
!      real specqqint(0:imx,hmax)
!      real speceeint(0:imx,hmax)
!      real specqcint(0:imx,hmax)
!      real spectpint(0:imx,hmax)
!      real spectsint(0:imx,hmax)
!      real specucint(0:imx,hmax)
!      real specfxint(0:imx,hmax)
!      real specepsint(0:imx,hmax)
!      real specepsrint(0:imx,hmax)
!      real specfyint(0:imx,hmax)
!      real specfzint(0:imx,hmax)
!      real spectsvint(0:imx,hmax)
!      real spectswint(0:imx,hmax)
!     
!      integer nfiles,rnk,rs
!      real stime
!
!      ucfluc=ufluc_s*cfluc
!      vcfluc=wfluc*cfluc
!      wcfluc=wfluc*cfluc
!      do i=1,imax
!       interp_s(i,:,:)=cfluc(i,:,:)
!      enddo
!      call inter_s_c_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
!      call deriv_c_s_m (imax,jmax*kmax/px,interp_c,deriv_s,dr)
!      do i=1,imax
!       dcfluc(i,:,:)=deriv_s(i,:,:)*mr_s(i)
!      enddo
!      do i=1,imax
!       interp_s(i,:,:)=ucfluc(i,:,:)
!      enddo
!      call inter_s_c_m(imax,jmax*kmax/px,interp_s,interp_c,dr)
!      call deriv_c_s_m (imax,jmax*kmax/px,interp_c,deriv_s,dr)
!      do i=1,imax
!       ducfluc(i,:,:)=deriv_s(i,:,:)*mr_s(i)
!      enddo      
!
!      call der1c_s_c6_m(imax,1,cm,dc,dr,Twb,Twt)
!      dc=dc*mr_c
!      call inter_c_s_m(imax,1,dc,dcc,dr)
!
!      if(rnk.eq.0) write(*,*) 'entering 2D spectra'
! 
!      stime = MPI_WTIME()
!      call spectra2Dfft(cfluc,cfluc,speccc,speccci) 
!      call spectra2Dfft(qfluc,qfluc,specqq,speccci) 
!      call spectra2Dfft(efluc,efluc,specee,speccci) 
!      call spectra2Dfft(qfluc,cfluc,specqc,speccci) 
!      call spectra2Dfft(wfluc,wfluc,specww,specwwi) 
!      call spectra2Dfft(ducfluc,cfluc,spectp,spectpi) 
!      call spectra2Dfft(vcfluc,cfluc,spectsv,spectsvi) 
!      call spectra2Dfft(wcfluc,cfluc,spectsw,spectswi) 
!      call spectra2Dfft(ufluc_s,cfluc,specuc,specuci) 
!      call spectra2Dfft(fxfluc,cfluc,specfxc,specfxci) 
!      call spectra2Dfft(fyfluc,cfluc,specfyc,specfyci) 
!      call spectra2Dfft(fzfluc,cfluc,specfzc,specfzci) 
!      call spectra2Dfft(dcfluc,dcfluc,speceps,specepsi) 
!      call spectra2Dfft(fxfluc,dcfluc,specepsr,specepsri) 
!       if(rnk.eq.0) write(*,*) 'Time in 2D spectra: ', MPI_WTIME()-stime
!
!      do i=0,imx
!       do j=1,jmax/2
!        do k=1,kmax/2
!         specfyci(i,j,k)=-specfyci(i,j,k)*(j-1)*2./1.5
!         spectsvi(i,j,k)=-spectsvi(i,j,k)*(j-1)*2./1.5
!         specfzci(i,j,k)=-specfzci(i,j,k)*(k-1)*2./4.
!         spectswi(i,j,k)=-spectswi(i,j,k)*(k-1)*2./4.
!        enddo
!       enddo
!      enddo
!
!      if(rnk.eq.0) write(*,*) 'entering averaging'
!
!      
!      call circint(speccc,specccint)
!      call circint(specqq,specqqint)
!      call circint(specee,speceeint)
!      call circint(specqc,specqcint)
!      call circint(specuc,specucint)
!      call circint(spectp,spectpint)
!      call circint(specfxc,specfxint)
!      call circint(speceps,specepsint)
!      call circint(specepsr,specepsrint)
!      call circint(specfyci,specfyint)
!      call circint(specfzci,specfzint)
!      call circint(spectsvi,spectsvint)
!      call circint(spectswi,spectswint)
! 
!      do i=1,imax
!       int_s(i,:)=specccint(i,:)
!      enddo
!      call deriv_s_c_m(imax,hmax,int_s,int_c,dr)
!      do h=1,hmax
!       int_c(:,h)=int_c(:,h)*mr_c(:)
!      enddo
!      call deriv_c_s_m(imax,hmax,int_c,der_s,dr)
!      do i=1,imax
!       int_s(i,:)=specfxint(i,:)
!      enddo
!      call inter_s_c_m(imax,hmax,int_s,int_c,dr)
!      call deriv_c_s_m (imax,hmax,int_c,der2_s,dr)
!      
!      if(rnk.eq.0) write(*,*) 'writing spectral budgets'
!
!      do i=1,imax
!       specbud(i,:,1)=specbud(i,:,1) - specucint(i,:)*dcc(i)/nfiles
!       specbud(i,:,2)=specbud(i,:,2) + (spectsvint(i,:)+spectswint(i,:))/nfiles
!       specbud(i,:,3)=specbud(i,:,3) - spectpint(i,:)/nfiles
!       do h=1,hmax
!        specbud(i,h,4)=specbud(i,h,4) - specccint(i,h)*(bp(h)**2.)/(Re*Pr*nfiles) 
!       enddo
!       specbud(i,:,5)=specbud(i,:,5) - specepsint(i,:)*2./(Re*Pr*nfiles)
!       specbud(i,:,6)=specbud(i,:,6) + der_s(i,:)*mr_s(i)*2./(Re*Pr*nfiles)
!       specbud(i,:,7)=specbud(i,:,7) - (specfyint(i,:)+specfzint(i,:))*2./(Re*Pr*Pl*nfiles) 
!       specbud(i,:,8)=specbud(i,:,8) + specepsrint(i,:)*2./(Re*Pr*Pl*nfiles) 
!       specbud(i,:,9)=specbud(i,:,9) - der2_s(i,:)*mr_s(i)*2./(Re*Pr*Pl*nfiles)
!      enddo
!
!      if(rnk.eq.0) write(*,*) 'calculating normalization factors'
!       
!      do rs=1,9
!       do i=1,imax
!        do h=1,hmax
!         specnorm(i,rs)=specnorm(i,rs)+(specbud(i,h,rs))*(bu(h)-bu(h-1))
!        enddo
!       enddo
!      enddo
!
!
!      return
!      end
!
!      subroutine circint(spectra,specint)
!      
!      implicit none
!      include 'par_post.txt'
!      include 'common.txt'
!
!      real spectra(0:imx,jmax/2,kmax/2)
!      real specint(0:imx,hmax)
!      integer nbin(0:imax,hmax)
!      real w1,w2,wr
!
!      specint=0.
!!      nbin=0      
!!      do i=0,imx
!!       do j=1,jmax/2
!!        do k=1,kmax/2
!!         w1=(j-1)*2./1.5
!!         w2=(k-1)*2./4.
!!         wr=sqrt(w1**2+w2**2)      
!!         do h=1,hmax
!!          if(wr.lt.bu(h).and.wr.gt.bu(h-1)) then
!!           nbin(i,h)=nbin(i,h)+1
!!           goto 133
!!          endif
!!         enddo
!!133      continue
!!        enddo
!!       enddo
!      do i=0,imx
!       do j=1,jmax/2
!        do k=1,kmax/2
!         w1=(j-1)*2./1.5
!         w2=(k-1)*2./4.
!         wr=sqrt(w1**2+w2**2)      
!         do h=1,hmax
!          if(wr.lt.bu(h).and.wr.ge.bu(h-1)) then
!           specint(i,h)=specint(i,h)+spectra(i,j,k)*bp(h)
!           goto 134
!          endif
!         enddo
!134      continue
!        enddo
!       enddo
!      enddo
!
!      return
!      end
!
!
!      subroutine spectra2Dfft(p1x,p2x,spectra,spectrai)
!
!      use decomp_2d
!      implicit none
!      include 'par_post.txt'
!      include 'common.txt'
!      include 'mpif.h'
!      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
!      real p2x(0:imax+1,jmax/p_row,kmax/p_col)
!      real p1y(0:imx,jmax,kmax/p_col)
!      real p2y(0:imx,jmax,kmax/p_col)
!      real p1c(0:imx,jmax,kmax)
!      real p2c(0:imx,jmax,kmax)
!      real send(0:imx,jmax,kmax)
!      complex four1(jmax,kmax)
!      complex four2(jmax,kmax)
!      complex spectrac1(0:imx,jmax,kmax)
!      complex spectrac2(0:imx,jmax,kmax)
!      real spectra2(0:imx,jmax,kmax)
!      real spectra2i(0:imx,jmax,kmax)
!      real spectra(0:imx,jmax/2,kmax/2)
!      real spectrai(0:imx,jmax/2,kmax/2)
!
!      integer lensav,lenwrk,ldim,seed
!      real, allocatable, dimension(:) :: work
!      real, allocatable, dimension(:) :: wsave
!      integer ier
!
!      p1c=0.
!      p2c=0.
!      call transpose_x_to_y(p1x,p1y)
!      call transpose_x_to_y(p2x,p2y)
!      do i=0,imx
!       do j=1,jmax
!        do k=1,kmax/p_col
!         kr=ystart(3)+k-1
!         p1c(i,j,kr)=p1y(i,j,k)
!         p2c(i,j,kr)=p2y(i,j,k)
!        enddo
!       enddo
!      enddo
!      call mpi_allreduce(p1c,send,(imx+1)*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
!      p1c=send
!      call mpi_allreduce(p2c,send,(imx+1)*jmax*kmax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
!      p2c=send
!
!      do i=0,imx
!
!       lensav = 2*(jmax+kmax)+int(log(real(jmax,kind=8)))+int(log(real(kmax,kind=8)))+8
!       lenwrk = 2*jmax*kmax
!
!       allocate(wsave(1:lensav))
!       allocate(work(1:lenwrk))
!
!       call zfft2i(jmax,kmax,wsave,lensav,ier)
!
!       do k=1,kmax
!        do j=1,jmax
!         four1(j,k) = cmplx(real(p1c(i,j,k),kind=8),0.0)
!         four2(j,k) = cmplx(real(p2c(i,j,k),kind=8),0.0)
!        enddo
!       enddo
!
!       call zfft2f(jmax,jmax,kmax,four1,wsave,lensav,work,lenwrk,ier)
!       call zfft2f(jmax,jmax,kmax,four2,wsave,lensav,work,lenwrk,ier)
!
!       deallocate(work)
!       deallocate(wsave)
!
!       spectrac1(i,:,:) = four1(:,:)
!       do j=1,jmax
!        do k=1,kmax
!         spectrac2(i,j,k) = conjg(four2(j,k))
!         spectra2(i,j,k) = real(spectrac1(i,j,k)*spectrac2(i,j,k))
!         spectra2i(i,j,k) = real(aimag(spectrac1(i,j,k)*spectrac2(i,j,k)))
!        enddo
!       enddo
!      enddo
!      do i=0,imx
!       do j=1,jmax/2
!        do k=1,kmax/2
!         spectra(i,j,k) = spectra2(i,j,k)
!         spectrai(i,j,k) = spectra2i(i,j,k)
!        enddo
!       enddo
!      enddo
!
!      return
!      end
!
!
!      subroutine mk_bin
!
!      implicit none
!      include 'par_post.txt'
!      include 'common.txt'
!
!      real bmin,bmax,db
!      real btemp(hmax)
!
!      bmin=1.
!      bu(0)=0.
!      db=2
!      do h=1,hmax
!!       btemp(h)=(db**h)*bmin
!       bu(h)=bu(h-1)+db
!      enddo
!      do h=1,hmax
!       bp(h)=(bu(h-1)+bu(h))/2.
!      enddo
! 
!      open(unit=1,file='bins')
!      do h=1,hmax
!       write(1,*) h,btemp(h),bu(h-1),bp(h)
!      enddo
!      close(1) 
! 
!      return
!      end
