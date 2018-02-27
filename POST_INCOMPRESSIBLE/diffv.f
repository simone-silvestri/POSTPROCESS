         enddo
             if (even.eq.'y') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j)/2+1,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
             if (even.eq.'n') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j+1)/2,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
          call dgetrs   (  'N',imax,2,mat,imax,indx,rhs2,imax,info)
          do i=1,imax
	  p1x(i,j  ,k)=rhs2(i,1)
	  p1x(i,j  ,k+1)=rhs2(i,2)
          enddo
        enddo
      enddo
        do k=kmax/p_col,kmax/p_col
         do  j =1,jmax/p_row
         do i=1,imax
	 rhs(i)   = p1x(i,j,k)
         enddo
             if (even.eq.'y') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j)/2+1,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
             if (even.eq.'n') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j+1)/2,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
           endif
          call dgetrs   (  'N',imax,1,mat,imax,indx,rhs,imax,info)
          do i=1,imax
	  p1x(i,j  ,k)=rhs(i)
          enddo
        enddo
       enddo




c********************************************************************
      call transpose_x_to_y(p1x,p1y)
      do k=1,kmax/p_col
        do i=0,imx
          do j=1,jmax
           ft(j)=p1y(i,j,k)
           enddo
           call vrfftb(1,jmax,ft,dft,1,wj)
           do j=1,jmax
         p1y(i,j,k)=ft(j)
           enddo
         enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
        do j=1,jmax/p_col
          do k=1,kmax
          fz(k)=p1z(i,j,k)
           enddo
           call vrfftb(1,kmax,fz,dfz,1,wk)
          do k=1,kmax
          p1z(i,j,k)=fz(k)
           enddo
         enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
            p(i,j,k)=p1x(i,j,k)
          enddo
        enddo
      enddo
 
      end











      SUBROUTINE DDGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
*     ..
*
*  Purpose
*  =======
*
*  DGTSV  solves the equation
*
*     A*X = B,
*
