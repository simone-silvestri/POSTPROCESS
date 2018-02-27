
      subroutine taubar(utmp,vtmp,wtmp,mu,srr,szz,stt,srz,srt,szt)
      use decomp_2d
      implicit none
      include 'par_post.txt'
      include 'common.txt'
      include 'mpif.h'   
 
      real d11,d12,d13,d21,d22,d23,d31,d32,d33,divtmp
      real d11p,d12p,d13p,d21p,d22p,d23p,d31p,d32p,d33p
      real S(3,3)
      real mumean(0:imax+1,jmax/p_row,kmax/p_col)
      real utmp(0:imax+1,jmax/p_row,kmax/p_col)
      real vtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real wtmp(0:imax+1,jmax/p_row,kmax/p_col)
      real hiu_s(imax,jmax/p_row,kmax/p_col),hiu_c(0:imax,jmax/p_row,kmax/p_col)
      real hiv_s(imax,jmax/p_row,kmax/p_col),hiv_c(0:imax,jmax/p_row,kmax/p_col)
      real hiw_s(imax,jmax/p_row,kmax/p_col),hiw_c(0:imax,jmax/p_row,kmax/p_col)
      real u_s(0:imax+1,jmax/p_row,kmax/p_col),u_c(0:imax+1,jmax/p_row,kmax/p_col)
      real v_s(0:imax+1,jmax/p_row,kmax/p_col),v_c(0:imax+1,jmax/p_row,kmax/p_col)
      real w_s(0:imax+1,jmax/p_row,kmax/p_col),w_c(0:imax+1,jmax/p_row,kmax/p_col)
      real dhiu_s(imax,jmax/p_row,kmax/p_col),dhiu_c(0:imax,jmax/p_row,kmax/p_col)
      real dhiv_s(imax,jmax/p_row,kmax/p_col),dhiv_c(0:imax,jmax/p_row,kmax/p_col)
      real dhiw_s(imax,jmax/p_row,kmax/p_col),dhiw_c(0:imax,jmax/p_row,kmax/p_col)
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax)
      real ft2(jmax),dft2(jmax),fz2(kmax),dfz2(kmax)
      real ft3(jmax),dft3(jmax),fz3(kmax),dfz3(kmax)
      real tmp (0:imx,jmax,kmax/p_col),tmpf (0:imx,jmax,kmax/p_col),tmpp(0:imx,jmax,kmax/p_col)
      real tmp2(0:imx,jmax,kmax/p_col),tmp2f(0:imx,jmax,kmax/p_col),tmp2p(0:imx,jmax,kmax/p_col)
      real tmp3(0:imx,jmax/p_col,kmax),tmp3f(0:imx,jmax/p_col,kmax),tmp3p(0:imx,jmax/p_col,kmax)
      real du_dt(0:imax+1,jmax/p_row,kmax/p_col),du2_dt(0:imax+1,jmax/p_row,kmax/p_col),dup_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dt(0:imax+1,jmax/p_row,kmax/p_col),dv2_dt(0:imax+1,jmax/p_row,kmax/p_col),dvp_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dt(0:imax+1,jmax/p_row,kmax/p_col),dw2_dt(0:imax+1,jmax/p_row,kmax/p_col),dwp_dt(0:imax+1,jmax/p_row,kmax/p_col)
      real du_dz(0:imax+1,jmax/p_row,kmax/p_col),du2_dz(0:imax+1,jmax/p_row,kmax/p_col),dup_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dv_dz(0:imax+1,jmax/p_row,kmax/p_col),dv2_dz(0:imax+1,jmax/p_row,kmax/p_col),dvp_dz(0:imax+1,jmax/p_row,kmax/p_col)
      real dw_dz(0:imax+1,jmax/p_row,kmax/p_col),dw2_dz(0:imax+1,jmax/p_row,kmax/p_col),dwp_dz(0:imax+1,jmax/p_row,kmax/p_col)
      integer rs


      real mu(0:imax+1,jmax/p_row,kmax/p_col)

      real srr(0:imax+1,jmax/p_row,kmax/p_col)
      real srt(0:imax+1,jmax/p_row,kmax/p_col)
      real srz(0:imax+1,jmax/p_row,kmax/p_col)
      real str(0:imax+1,jmax/p_row,kmax/p_col)
      real stt(0:imax+1,jmax/p_row,kmax/p_col)
      real stz(0:imax+1,jmax/p_row,kmax/p_col)
      real szr(0:imax+1,jmax/p_row,kmax/p_col)
      real szt(0:imax+1,jmax/p_row,kmax/p_col)
      real szz(0:imax+1,jmax/p_row,kmax/p_col)
      real term1(0:imax+1,jmax/p_row,kmax/p_col)
  
      

 
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         hiv_s(i,j,k) = vtmp(i,j,k)     
         hiw_s(i,j,k) = wtmp(i,j,k)    
        enddo
        do i=0,imax
         hiu_c(i,j,k) = utmp(i,j,k)    
        enddo
       enddo
      enddo


      call inter_c_s_m(imax,jmax*kmax/px,hiu_c,hiu_s,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiv_s,hiv_c,dr)
      call inter_s_c_m(imax,jmax*kmax/px,hiw_s,hiw_c,dr)



      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         u_s(i,j,k) = hiu_s(i,j,k)
         v_s(i,j,k) = hiv_s(i,j,k)
         w_s(i,j,k) = hiw_s(i,j,k)
        enddo
        do i=0,imax
         u_c(i,j,k)  = hiu_c(i,j,k)
         v_c(i,j,k)  = hiv_c(i,j,k)
         w_c(i,j,k)  = hiw_c(i,j,k)
        enddo
       enddo
      enddo


       call deriv_c_s_m(imax,jmax*kmax/px,hiu_c,dhiu_s,dr)
       call deriv_c_s_m(imax,jmax*kmax/px,hiv_c,dhiv_s,dr)
       call deriv_c_s_m(imax,jmax*kmax/px,hiw_c,dhiw_s,dr)

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
 
         divtmp= (d11+d22+d33)*2./3.     

          srr(i,j,k) =mu(i,j,k)*(2* d11-divtmp)
          stt(i,j,k) =mu(i,j,k)*(2* d22-divtmp)
          szz(i,j,k) =mu(i,j,k)*(2* d33-divtmp)
          srt(i,j,k) =mu(i,j,k)*(d12 + d21)
          srz(i,j,k) =mu(i,j,k)*(d13 + d31)
          str(i,j,k) =mu(i,j,k)*(d21 + d12)
          stz(i,j,k) =mu(i,j,k)*(d32 + d23)
          szr(i,j,k) =mu(i,j,k)*(d31 + d13)
          szt(i,j,k) =mu(i,j,k)*(d32 + d23)
    
         enddo
       enddo
      enddo



      end

      subroutine bound(u,v,w,e,c,rank)
      include 'par_post.txt'
      include 'common.txt'
      real u(0:i1,jmax/p_row,kmax/p_col)
      real v(0:i1,jmax/p_row,kmax/p_col)
      real w(0:i1,jmax/p_row,kmax/p_col)
      real c(0:i1,jmax/p_row,kmax/p_col)
      real e(0:i1,jmax/p_row,kmax/p_col)
      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  u(imax,j,k)=0.
	  v(i1  ,j,k)=-v(imax,j,k)	
	  w(i1  ,j,k)=-w(imax,j,k)	
	  e(i1,  j,k)=-e(imax,j,k)
	  c(i1,  j,k)=-c(imax,j,k)
	  u(0   ,j,k)= 0!-u(1,j,k)
	  v(0   ,j,k)= -v(1,j,k)
	  w(0   ,j,k)= -w(1,j,k)
	  e(0   ,j,k)= (e(1,j,k)-1.0)/(rp(1)-0.5)*(rp(0) - rp(1)) + e(1,j,k)
	  c(0   ,j,k)= (c(1,j,k)-1.0)/(rp(1)-0.5)*(rp(0) - rp(1)) + c(1,j,k)
	enddo
       enddo
      end

      subroutine spectrum(nmax,x,spec)
      implicit none
      integer nmax,istap
      real     x(nmax),spec(nmax/2)
      complex  a(nmax),ii,work(2*nmax)
      complex  a2(nmax)

c   Declarations  for FFT
      integer  ifax(13),inc,jump,lot,n,i,isign
      real     trigs(2*nmax),t,pi

      ii = ( 0.,1.)
      call cftfax(nmax,ifax,trigs)
      if (ifax(1).eq. -99) stop 'No factor 2**p * 3**q * 5**r'

      pi =4.*atan(1.)

c
c     read data
c

      do i=1,nmax
	a(i) = x(i)
      enddo

c isign = -1   naar Fourier Ruimte
c isign =  1   Terug

      isign = -1
      inc   =  1
      jump  =  1
      lot   =  1
      n     =  nmax

      call cfft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
      
      do i=1,nmax
      a2(i) = (a(i)*conjg(a(i)))/(nmax**2)
      enddo 
      
c
c     Power spectrum
c
      do i=1,nmax/2
      spec(i)= real(a2(i))
      enddo
      end
      
c******************************************************************************

      SUBROUTINE CFFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC
C              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.
C              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE
C              PACKAGE RETURNS A SET OF COMPLEX FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS
C              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.
C
C              THE PACKAGE CFFT99 CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE CFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO CFFT99
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINE CFFT99
C                THE ACTUAL TRANSFORM ROUTINE ROUTINE, CABABLE OF
C                PERFORMING BOTH THE TRANSFORM AND ITS INVERSE.
C                HOWEVER, AS THE TRANSFORMS ARE NOT NORMALIZED,
C                THE APPLICATION OF A TRANSFORM FOLLOWED BY ITS
C                INVERSE WILL YIELD THE ORIGINAL VALUES MULTIPLIED
C                BY N.
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF
C              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER
C              COEFFICIENT VECTORS OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(2*N)
C                   COMPLEX A(...), WORK(...)
C
C                   CALL CFTFAX (N, IFAX, TRIGS)
C                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND
C              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR
C              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED
C              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.
C              SEE THE INDIVIDUAL WRITE-UPS FOR CFTFAX AND
C              CFFT99 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.  IT WAS
C              FURTHER MODIFIED FOR THE FULLY COMPLEX CASE BY DAVE
C              FULKER IN NOVEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,
C              PROVIDED THAT N IS NOT CHANGED.
C
C ARGUMENT     IFAX(13),TRIGS(2*N)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF
C               THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.
C
C              TRIGS
C               A REAL ARRAY OF DIMENSION 2*N
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS
C               GREATER THAN 5, IFAX(1) IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE CFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS (UNNORMALIZED) COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS.  GIVEN A SET OF COMPLEX GRIDPOINT
C              VECTORS, THE PACKAGE RETURNS A SET OF
C              COMPLEX FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE A
C              NUMBER HAVING NO PRIME FACTORS OTHER THAN
C              2, 3, AND 5.  THIS ROUTINE IS
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     COMPLEX A(N*INC+(LOT-1)*JUMP), WORK(N*LOT)
C DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13)
C
C ARGUMENTS
C
C ON INPUT     A
C               A COMPLEX ARRAY OF LENGTH N*INC+(LOT-1)*JUMP CONTAINING
C               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY
C               OVERWRITTEN BY THE RESULTS.
C
C               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF
C               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS
C               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT
C               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER"S"
C               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS
C               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE
C               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A
C               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF
C               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE
C               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND
C               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND
C               IMAGINARY PORTIONS EXPLICITLY.
C
C              WORK
C               A COMPLEX WORK ARRAY OF LENGTH N*LOT OR A REAL ARRAY
C               OF LENGTH 2*N*LOT.  SEE N.B. ABOVE.
C
C              TRIGS
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C
C               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED
C               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED
C               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.
C
C              INC
C               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENT
C               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR
C               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS
C               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              LOT
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C
C ON OUTPUT    A
C               IF ISIGN = -1, AND LOT GRIDPOINT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N COEFFICIENT VALUES:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 C(K) = SUM(J=0,...,N-1)( G(J)*EXP(-2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               IF ISIGN = +1, AND LOT COEFFICIENT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N GRIDPOINT VALUES:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 G(J) = SUM(K=0,...,N-1)( G(K)*EXP(+2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               A CALL WITH ISIGN=-1 FOLLOWED BY A CALL WITH ISIGN=+1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED
C               BY THE FACTOR N.
C
C
C EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN
C               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL
C               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,
C               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE
C               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE
C               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE
C               ORTHOGONAL DIRECTION.
C
C               COMPLEX A(66,9), WORK(64,9)
C               REAL TRIGS1(128), TRIGS2(18)
C               INTEGER IFAX1(13), IFAX2(13)
C
C               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:
C
C               CALL CFTFAX(64, IFAX1, TRIGS1)
C               CALL CFTFAX( 9, IFAX2, TRIGS2)
C
C               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE
C               STORED IN MEMORY AS FOLLOWS (USING U AND V TO
C               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND
C               ASSUMING CONVENTIONAL FORTRAN STORAGE):
C
C   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS,
C
C   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS,
C
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C
C   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS.
C
C               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE
C               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 9, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST
C               DIMENSION OF A).
C
C               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)
C
C               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT
C               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFIN
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 64, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 1
C
C               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)
C
C               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONA
C               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT
C               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN
C               WITH ISIGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF
C               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).
C
C
C-----------------------------------------------------------------------
      DIMENSION A(1),WORK(1),TRIGS(1),IFAX(1)
C
C     SUBROUTINE "CFFT99" - MULTIPLE FAST COMPLEX FOURIER TRANSFORM
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE N*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL.
C
C
C
C
C
      NN = N+N
      INK=INC+INC
      JUM = JUMP+JUMP
      NFAX=IFAX(1)
      JNK = 2
      JST = 2
      IF (ISIGN.GE.0) GO TO 30
C
C     THE INNERMOST TEMPERTON ROUTINES HAVE NO FACILITY FOR THE
C     FORWARD (ISIGN = -1) TRANSFORM.  THEREFORE, THE INPUT MUST BE
C     REARRANGED AS FOLLOWS:
C
C     THE ORDER OF EACH INPUT VECTOR,
C
C     G(0), G(1), G(2), ... , G(N-2), G(N-1)
C
C     IS REVERSED (EXCLUDING G(0)) TO YIELD
C
C     G(0), G(N-1), G(N-2), ... , G(2), G(1).
C
C     WITHIN THE TRANSFORM, THE CORRESPONDING EXPONENTIAL MULTIPLIER
C     IS THEN PRECISELY THE CONJUGATE OF THAT FOR THE NORMAL
C     ORDERING.  THUS THE FORWARD (ISIGN = -1) TRANSFORM IS
C     ACCOMPLISHED
C
C     FOR NFAX ODD, THE INPUT MUST BE TRANSFERRED TO THE WORK ARRAY,
C     AND THE REARRANGEMENT CAN BE DONE DURING THE MOVE.
C
      JNK = -2
      JST = NN-2
      IF (MOD(NFAX,2).EQ.1) GOTO 40
C
C     FOR NFAX EVEN, THE REARRANGEMENT MUST BE APPLIED DIRECTLY TO
C     THE INPUT ARRAY.  THIS CAN BE DONE BY SWAPPING ELEMENTS.
C
      IBASE = 1
      ILAST = (N-1)*INK
      NH = N/2
      DO 20 L=1,LOT
      I1 = IBASE+INK
      I2 = IBASE+ILAST
CDIR$ IVDEP
      DO 10 M=1,NH
C     SWAP REAL AND IMAGINARY PORTIONS
      HREAL = A(I1)
      HIMAG = A(I1+1)
      A(I1) = A(I2)
      A(I1+1) = A(I2+1)
      A(I2) = HREAL
      A(I2+1) = HIMAG
      I1 = I1+INK
      I2 = I2-INK
   10 CONTINUE
      IBASE = IBASE+JUM
   20 CONTINUE
      GOTO 100
C
   30 CONTINUE
      IF (MOD(NFAX,2).EQ.0) GOTO 100
C
   40 CONTINUE
C
C     DURING THE TRANSFORM PROCESS, NFAX STEPS ARE TAKEN, AND THE
C     RESULTS ARE STORED ALTERNATELY IN WORK AND IN A.  IF NFAX IS
C     ODD, THE INPUT DATA ARE FIRST MOVED TO WORK SO THAT THE FINAL
C     RESULT (AFTER NFAX STEPS) IS STORED IN ARRAY A.
C
      IBASE=1
      JBASE=1
      DO 60 L=1,LOT
C     MOVE REAL AND IMAGINARY PORTIONS OF ELEMENT ZERO
      WORK(JBASE) = A(IBASE)
      WORK(JBASE+1) = A(IBASE+1)
      I=IBASE+INK
      J=JBASE+JST
CDIR$ IVDEP
      DO 50 M=2,N
C     MOVE REAL AND IMAGINARY PORTIONS OF OTHER ELEMENTS (POSSIBLY IN
C     REVERSE ORDER, DEPENDING ON JST AND JNK)
      WORK(J) = A(I)
      WORK(J+1) = A(I+1)
      I=I+INK
      J=J+JNK
   50 CONTINUE
      IBASE=IBASE+JUM
      JBASE=JBASE+NN
   60 CONTINUE
C
  100 CONTINUE
C
C     PERFORM THE TRANSFORM PASSES, ONE PASS FOR EACH FACTOR.  DURING
C     EACH PASS THE DATA ARE MOVED FROM A TO WORK OR FROM WORK TO A.
C
C     FOR NFAX EVEN, THE FIRST PASS MOVES FROM A TO WORK
      IGO = 110
C     FOR NFAX ODD, THE FIRST PASS MOVES FROM WORK TO A
      IF (MOD(NFAX,2).EQ.1) IGO = 120
      LA=1
      DO 140 K=1,NFAX
      IF (IGO.EQ.120) GO TO 120
  110 CONTINUE
      CALL VPASSM(A(1),A(2),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUM,NN,LOT,N,IFAX(K+1),LA)
      IGO=120
      GO TO 130
  120 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(1),A(2),TRIGS,
     *    2,INK,NN,JUM,LOT,N,IFAX(K+1),LA)
      IGO=110
  130 CONTINUE
      LA=LA*IFAX(K+1)
  140 CONTINUE
C
C     AT THIS POINT THE FINAL TRANSFORM RESULT IS STORED IN A.
C
      RETURN
      END
      SUBROUTINE CFTFAX(N,IFAX,TRIGS)
      DIMENSION IFAX(13),TRIGS(1)
C
C     THIS ROUTINE WAS MODIFIED FROM TEMPERTON"S" ORIGINAL
C     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING
C     ORDER, AND THERE ARE NONE OF THE ORIGINAL 'MODE' OPTIONS.
C
C ON INPUT     N
C               THE LENGTH OF EACH COMPLEX TRANSFORM TO BE PERFORMED
C
C               N MUST BE GREATER THAN 1 AND CONTAIN NO PRIME
C               FACTORS GREATER THAN 5.
C
C ON OUTPUT    IFAX
C               IFAX(1)
C                 THE NUMBER OF FACTORS CHOSEN OR -99 IN CASE OF ERROR
C               IFAX(2) THRU IFAX( IFAX(1)+1 )
C                 THE FACTORS OF N IN THE FOLLOWIN ORDER:  APPEARING
C                 FIRST ARE AS MANY FACTORS OF 4 AS CAN BE OBTAINED.
C                 SUBSEQUENT FACTORS ARE PRIMES, AND APPEAR IN
C                 ASCENDING ORDER, EXCEPT FOR MULTIPLE FACTORS.
C
C              TRIGS
C               2N SIN AND COS VALUES FOR USE BY THE TRANSFORM ROUTINE
C
      CALL FACT(N,IFAX)
      K = IFAX(1)
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 ) THEN
        WRITE(6,*) ' FFTFAX -- INVALID N'
      ENDIF
      CALL CFTRIG (N, TRIGS)
      RETURN
      END
      SUBROUTINE FACT(N,IFAX)
C     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4
      DIMENSION IFAX(13)
      IF (N.GT.1) GO TO 10
      IFAX(1) = 0
      IF (N.LT.1) IFAX(1) = -99
      RETURN
   10 NN=N
      K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      MAX = SQRT(FLOAT(NN))
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 IF (L.GT.MAX) GO TO 75
      L=L+INC
      INC=6-INC
      GO TO 60
   75 K = K+1
      IFAX(K) = NN
   80 IFAX(1)=K-1
C     IFAX(1) NOW CONTAINS NUMBER OF FACTORS
      RETURN
      END
      SUBROUTINE CFTRIG(N,TRIGS)
      DIMENSION TRIGS(1)
      PI=2.0*ASIN(1.0)
      DEL=(PI+PI)/FLOAT(N)
      L=N+N
      DO 10 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
C
C     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
C     PERFORMS ONE PASS THROUGH DATA
C     AS PART OF MULTIPLE COMPLEX (INVERSE) FFT ROUTINE
C     A IS FIRST REAL INPUT VECTOR
C     B IS FIRST IMAGINARY INPUT VECTOR
C     C IS FIRST REAL OUTPUT VECTOR
C     D IS FIRST IMAGINARY OUTPUT VECTOR
C     TRIGS IS PRECALCULATED TABLE OF SINES & COSINES
C     INC1 IS ADDRESSING INCREMENT FOR A AND B
C     INC2 IS ADDRESSING INCREMENT FOR C AND D
C     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
C     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
C     LOT IS THE NUMBER OF VECTORS
C     N IS LENGTH OF VECTORS
C     IFAC IS CURRENT FACTOR OF N
C     LA IS PRODUCT OF PREVIOUS FACTORS
C
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,
     *     SIN60/0.866025403784437/
C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
C     CODING FOR FACTOR 2
C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 3
C
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     *    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 4
C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 5
C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END

      subroutine four1(n,a,da,dx)
      implicit none
      integer n,i, j,jwav
      real a(n),da(n),dx,h,rn
      real wk(2*n+15)
      call vrffti(n,wk)
c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo
      call vrfftb(1,n,da,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
      enddo
      end



      subroutine four12(n,a,da,dda,dx)
      implicit none
      integer n,i, j,jwav
      real a(n),da(n),dx,h,rn,dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)
c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      dda = da
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo
      do j=1,n
         jwav = j/2
	 dda(j)=-dda(j)*(jwav)**2
      enddo
      call vrfftb(1,n,da,a,1,wk)
      call vrfftb(1,n,dda,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
	dda(i)=dda(i)*rn*rn
      enddo
      end
      subroutine four1r(n,itr,a,da,dx)
      implicit none
      integer n,i, j,jwav,itr
      real a(n),da(n),dx,h,rn
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo

      do j=itr,n
      da(j)=0
      enddo

      call vrfftb(1,n,da,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
      enddo
      end

      subroutine four2(n,a,dda,dx)
      integer n
      real a(n),dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      do i=1,n
	 dda(i)=a(i)
      enddo
      call vrfftf(1,n,dda,a,1,wk)

      do j=1,n
         jwav = j/2  
         dda(j) =- dda(j)*jwav**2
      enddo

c
c    Back to physical space
c
      call vrfftb(1,n,dda,a,1,wk)
      rn = (8*atan(1.)/(n*dx))**2

      do j=1,n
	   dda(j)=dda(j)*rn
      enddo
      end

      subroutine four2r(n,itr,a,dda,dx)
      integer n,itr
      real a(n),dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      dda =  0
      do i=1,n
	 dda(i)=a(i)
      enddo
      call vrfftf(1,n,dda,a,1,wk)

      do j=1,n
         jwav = j/2  
         dda(j) =- dda(j)*jwav**2
      enddo
      do j=itr,n
       dda(j)=0
      enddo

c
c    Back to physical space
c
      call vrfftb(1,n,dda,a,1,wk)
c    
      rn = (8*atan(1.)/(n*dx))**2

      do j=1,n
	   dda(j)=dda(j)*rn
      enddo
c
c   all done
c 
      end

       function pimach()
       pimach = 4.*atan(1.)
       end


      SUBROUTINE VRADB2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,2,L1)    ,CH(MDIMC,IDO,L1,2),
     1                WA1(IDO)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)-CC(M,IDO,2,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+CC(M,IC-1,2,K)
            CH(M,I,K,1) = CC(M,I,1,K)-CC(M,IC,2,K)
            CH(M,I-1,K,2) = WA1(I-2)*(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
     1      -WA1(I-1)*(CC(M,I,1,K)+CC(M,IC,2,K))
            CH(M,I,K,2) = WA1(I-2)*(CC(M,I,1,K)+CC(M,IC,2,K))+WA1(I-1)
     1      *(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
          DO 1003 M=1,MP
         CH(M,IDO,K,1) = CC(M,IDO,1,K)+CC(M,IDO,1,K)
         CH(M,IDO,K,2) = -(CC(M,1,2,K)+CC(M,1,2,K))
 1003     CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADB3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,3,L1)    ,CH(MDIMC,IDO,L1,3),
     1                WA1(IDO)   ,WA2(IDO)
      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   -(2.*TAUI)*CC(M,1,3,K)
         CH(M,1,K,3) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   +2.*TAUI*CC(M,1,3,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2) = WA1(I-2)*
     1 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     * (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     2                   -WA1(I-1)*
     3 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     * (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,2) = WA1(I-2)*
     4 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     5                  +WA1(I-1)*
     6 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
              CH(M,I-1,K,3) = WA2(I-2)*
     7 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     8                      -WA2(I-1)*
     9 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,3) = WA2(I-2)*
     1 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     2                 +WA2(I-1)*
     3 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADB4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,4,L1)  ,CH(MDIMC,IDO,L1,4)    ,
     1                WA1(IDO)  ,WA2(IDO)  ,WA3(IDO)
      SQRT2=SQRT(2.)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,3) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   -(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,1) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   +(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,4) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   +(CC(M,1,3,K)+CC(M,1,3,K))
         CH(M,1,K,2) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   -(CC(M,1,3,K)+CC(M,1,3,K))
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = (CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      +(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = (CC(M,I,1,K)-CC(M,IC,4,K))
     1      +(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2)=WA1(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      -(CC(M,I,3,K)+CC(M,IC,2,K)))-WA1(I-1)
     1      *((CC(M,I,1,K)+CC(M,IC,4,K))+(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,2)=WA1(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      +(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA1(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))-(CC(M,I,3,K)+CC(M,IC,2,K)))
            CH(M,I-1,K,3)=WA2(I-2)*((CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      -(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-WA2(I-1)
     1      *((CC(M,I,1,K)-CC(M,IC,4,K))-(CC(M,I,3,K)-CC(M,IC,2,K)))
            CH(M,I,K,3)=WA2(I-2)*((CC(M,I,1,K)-CC(M,IC,4,K))
     1      -(CC(M,I,3,K)-CC(M,IC,2,K)))+WA2(I-1)
     1      *((CC(M,I-1,1,K)+CC(M,IC-1,4,K))-(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K)))
            CH(M,I-1,K,4)=WA3(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      +(CC(M,I,3,K)+CC(M,IC,2,K)))-WA3(I-1)
     1     *((CC(M,I,1,K)+CC(M,IC,4,K))-(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,4)=WA3(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      -(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA3(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))+(CC(M,I,3,K)+CC(M,IC,2,K)))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
               DO 1003 M=1,MP
         CH(M,IDO,K,1) = (CC(M,IDO,1,K)+CC(M,IDO,3,K))
     1   +(CC(M,IDO,1,K)+CC(M,IDO,3,K))
         CH(M,IDO,K,2) = SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   -(CC(M,1,2,K)+CC(M,1,4,K)))
         CH(M,IDO,K,3) = (CC(M,1,4,K)-CC(M,1,2,K))
     1   +(CC(M,1,4,K)-CC(M,1,2,K))
         CH(M,IDO,K,4) = -SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   +(CC(M,1,2,K)+CC(M,1,4,K)))
 1003          CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADB5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,5,L1)    ,CH(MDIMC,IDO,L1,5),
     1             WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
      DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)+2.*CC(M,IDO,4,K)
         CH(M,1,K,2) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))-(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
         CH(M,1,K,3) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))-(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,4) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))+(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,5) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))+(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
 1001          CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
      DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +(CC(M,I-1,5,K)+CC(M,IC-1,4,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +(CC(M,I,5,K)-CC(M,IC,4,K))
            CH(M,I-1,K,2) = WA1(I-2)*((CC(M,I-1,1,K)+TR11*
     1      (CC(M,I-1,3,K)+CC(M,IC-1,2,K))+TR12
     1      *(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA1(I-1)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))+(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,2) = WA1(I-2)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)
     1      -CC(M,IC,2,K))+TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI11*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))+TI12
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))+WA1(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K))+TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))
     1      -(TI11*(CC(M,I,3,K)+CC(M,IC,2,K))+TI12
     1      *(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,3) = WA2(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1     -WA2(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,3) = WA2(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA2(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,4) = WA3(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA3(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,4) = WA3(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA3(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,5) = WA4(I-2)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA4(I-1)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,5) = WA4(I-2)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA4(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADBG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
     *                 MDIMC,WA)
      DIMENSION    CH(MDIMC,IDO,L1,IP)    ,CC(MDIMC,IDO,IP,L1) ,
     1           C1(MDIMC,IDO,L1,IP)     ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)       ,WA(IDO)
      TPI=2.*PIMACH(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            DO 1001 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1001       CONTINUE
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            DO 1004 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1004       CONTINUE
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            DO 1007 M=1,MP
            CH(M,1,K,J) = CC(M,IDO,J2-2,K)+CC(M,IDO,J2-2,K)
            CH(M,1,K,JC) = CC(M,1,J2-1,K)+CC(M,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               DO 1009 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1009          CONTINUE
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               DO 1013 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1013          CONTINUE
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            DO 1017 M=1,MP
            C2(M,IK,L) = CH2(M,IK,1)+AR1*CH2(M,IK,2)
            C2(M,IK,LC) = AI1*CH2(M,IK,IP)
 1017       CONTINUE
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               DO 1018 M=1,MP
               C2(M,IK,L) = C2(M,IK,L)+AR2*CH2(M,IK,J)
               C2(M,IK,LC) = C2(M,IK,LC)+AI2*CH2(M,IK,JC)
 1018          CONTINUE
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            DO 1021 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+CH2(M,IK,J)
 1021       CONTINUE
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            DO 1023 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)-C1(M,1,K,JC)
            CH(M,1,K,JC) = C1(M,1,K,J)+C1(M,1,K,JC)
 1023       CONTINUE
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               DO 1025 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               DO 1029 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1029          CONTINUE
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         DO 1033 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1033    CONTINUE
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            DO 1034 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)
 1034       CONTINUE
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               DO 1036 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1036          CONTINUE
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1040 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1040          CONTINUE
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
      SUBROUTINE VRADF2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,2,L1)  ,CC(MDIMC,IDO,L1,2)     ,
     1                WA1(IDO)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))
            CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-CC(M,I,K,1)
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         DO 1006 M=1,MP
         CH(M,1,2,K) = -CC(M,IDO,K,2)
         CH(M,IDO,1,K) = CC(M,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
         CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR*
     1      (CC(M,1,K,2)+CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))
            CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADF4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION    CC(MDIMC,IDO,L1,4)   ,CH(MDIMC,IDO,4,L1)     ,
     1                WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)
      HSQT2=SQRT(2.)/2.
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4))
     1      +(CC(M,1,K,1)+CC(M,1,K,3))
         CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3))
     1      -(CC(M,1,K,2)+CC(M,1,K,4))
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
         CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+
     1       WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         DO 1006 M=1,MP
            CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,1)
            CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)-
     1       CC(M,IDO,K,4)))
            CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))-
     1       CC(M,IDO,K,3)
            CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,L1,5)    ,CH(MDIMC,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADFG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,MDIMC,WA)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION     CH(MDIMC,IDO,L1,IP)   ,CC(MDIMC,IDO,IP,L1)  ,
     1            C1(MDIMC,IDO,L1,IP)    ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)           ,WA(IDO)
      TPI=2.*PIMACH(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         DO 1001 M=1,MP
         CH2(M,IK,1) = C2(M,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            DO 1002 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               DO 1004 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1008 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               DO 1012 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               DO 1016 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         DO 1020 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            DO 1022 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
            C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            DO 1024 M=1,MP
            CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
            CH2(M,IK,LC) = AI1*C2(M,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               DO 1025 M=1,MP
               CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
               CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            DO 1028 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            DO 1030 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            DO 1033 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            DO 1036 M=1,MP
            CC(M,IDO,J2-2,K) = CH(M,1,K,J)
            CC(M,1,J2-1,K) = CH(M,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               DO 1038 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               DO 1042 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE VRFFTB(M,N,R,RT,MDIMR,WSAVE)
C***BEGIN PROLOGUE  VRFFTB
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Backward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTB computes the synthesis (backward transform) of a
C  number of real periodic sequences from their Fourier coefficients. 
C  Specifically, for each set of independent Fourier coefficients
C  F(K), the corresponding real periodic sequence is computed. 
C
C  The array WSAVE which is used by subroutine VRFFTB must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sets of coefficients.
C
C  N       the length of the sequences of coefficients to be 
C          transformed.  The method is most efficient when N is a
C          product of small primes, however n may be any positive 
C          integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          coefficients to be transformed.  Each set of coefficients
C          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
C          the I-th set of independent Fourier coefficients is stored
C
C                R(I,1) = REAL( F(I,0) ),
C
C                R(I,2*K) = REAL( F(I,K) )
C
C                R(I,2*K+1) = IMAG( F(I,K) )
C
C                   for K = 1, 2, . . . , M-1,
C
C                and, when N is even,
C
C                R(I,N) = REAL( F(I,N/2) ).
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTB and VRFFTB.
C
C  Output Parameters
C
C  R       contains M real periodic sequences corresponding to the given
C          coefficients.  Specifically, the I-th row of R contains the 
C          real periodic sequence corresponding to the I-th set of
C          independent Fourier coefficients F(I,K) stored as
C
C               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
C
C               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
C                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is even, and
C
C               X(I,J) = SQRT(1/N)* F(I,0) +
C                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is odd.
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTB is a straightforward extension of the subprogram RFFTB to
C  handle M simultaneous sequences.  RFFTB was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTB1
C***END PROLOGUE  VRFFTB
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION     R(MDIMR,N),RT(MDIMR,N),WSAVE(N+15)
      IF (N .EQ. 1) RETURN
      CALL VRFTB1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFFTF (M,N,R,RT,MDIMR,WSAVE)
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward 
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M 
C          input sequences.  Specifically, row I of R, R(I,J), 
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = REAL( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = REAL( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = REAL( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       R(MDIMR,N)  ,RT(MDIMR,N)    ,WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFFTI (N,WSAVE)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (N .EQ. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFTB1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N), C(MDIMC,N), WA(N) ,FAC(15)
      NF = FAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADB4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL VRADB4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL VRADB2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 105
  104    CALL VRADB2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL VRADB3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 108
  107    CALL VRADB3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
      CALL VRADB5 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110 CALL VRADB5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL VRADBG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         GO TO 114
  113    CALL VRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 0) GO TO 118
      DO 117 J=1,N
      DO 117 I=1,M
         C(I,J) = SCALE*CH(I,J)
  117 CONTINUE
      RETURN
  118 DO 119 J=1,N
      DO 119 I=1,M
         C(I,J)=SCALE*C(I,J)
  119 CONTINUE
      RETURN
      END
      SUBROUTINE VRFTF1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N) ,C(MDIMC,N)  ,WA(N)   ,FAC(15)
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 1) GO TO 113
      DO 112 J=1,N
      DO 112 I=1,M
         C(I,J) = SCALE*CH(I,J)
  112 CONTINUE
      RETURN
  113 DO 114 J=1,N
      DO 114 I=1,M
         C(I,J)=SCALE*C(I,J)
  114 CONTINUE
      RETURN
      END
      SUBROUTINE VRFTI1 (N,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WA(N)      ,FAC(15)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 2.*PIMACH(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

      SUBROUTINE TRID (MR,A,B,C,Y,D)
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,Y(1)       ,
     1                D(1)
      call dgtsv(mr,1,a(2),b,c,y,mr,info)
      if (info.ne.0) stop 'error in tridagional systems ' 
      END  


      SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
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
*  where A is an n by n tridiagonal matrix, by Gaussian elimination with
*  partial pivoting.
*
*  Note that the equation  A'*X = B  may be solved by interchanging the '
*  order of the arguments DU and DL.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, DL must contain the (n-1) sub-diagonal elements of
*          A.
*
*          On exit, DL is overwritten by the (n-2) elements of the
*          second super-diagonal of the upper triangular matrix U from
*          the LU factorization of A, in DL(1), ..., DL(n-2).
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, D must contain the diagonal elements of A.
*
*          On exit, D is overwritten by the n diagonal elements of U.
*
*  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, DU must contain the (n-1) super-diagonal elements
*          of A.
*
*          On exit, DU is overwritten by the (n-1) elements of the first
*          super-diagonal of U.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N by NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N by NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
*               has not been computed.  The factorization has not been
*               completed unless i = N.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   FACT, TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Subroutines ..
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( NRHS.EQ.1 ) THEN
         DO 10 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
*
*              No row interchange required
*
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
*
*              Interchange rows I and I+1
*
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
   10    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      ELSE
         DO 40 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
*
*              No row interchange required
*
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
*
*              Interchange rows I and I+1
*
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            END IF
   40    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      END IF
*
*     Back solve with the matrix U from the factorization.
*
      IF( NRHS.LE.2 ) THEN
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 )
     $      B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*
     $                  B( I+2, J ) ) / D( I )
   80    CONTINUE
         IF( J.LT.NRHS ) THEN
            J = J + 1
            GO TO 70
         END IF
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 )
     $         B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /
     $                       D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*
     $                     B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of DGTSV
*
      END
      subroutine inter_c_s(n,f,df,dx)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call interpol_c_s6('n',n,f,df,dx,t1,a1,r)
      end 

      subroutine inter_s_c(n,f,df,dx)
      real F(N),dF(0:N),dx,t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call interpol_s_c6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine deriv_c_s(n,f,df,dx)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call der1_c_s6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine deriv_s_c(n,f,df,dx)
      real F(N),dF(0:N),dx ,t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call der1_s_c6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine inter_s_c_m(n,nrhs,f,df,dx)
      real F(N,nrhs),dF(0:N,nrhs),dx,t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      write(6,*) 'dentro inter_s_c_m'
      call interpol_s_c6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      write(6,*) 'fuori interpol_s_c6_m'
      end

      subroutine inter_c_s_m(n,nrhs,f,df,dx)
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call interpol_c_s6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end 

      subroutine deriv_c_s_m(n,nrhs,f,df,dx)
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call der1_c_s6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end

      subroutine deriv_s_c_m(n,nrhs,f,df,dx)
      real F(N),dF(0:N),dx ,t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call der1_s_c6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end

      subroutine der2_s_v(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=0.5
       val = 0.0
       rhs(1)=(16.*val/5. -9./2.*f(1)+f(2)+(3./10.)*f(3))/(dx*dx)
       b(n)=1.
       c(n)=0.
       a(n)=0.5
       val = 1.0
       rhs(n)=(16.*val/5. -9./2.*f(n)+f(n-1)+(3./10.)*f(n-2))/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
      end

      subroutine der2_s_n(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=-11./23
       val = 0.0
       rhs(1)=-24.*val/(23.*dx) +( -36./23.*f(1)+48*f(2)/23.-12.*f(3)/23.)/(dx*dx)

       b(n)=1.
       c(n)=0.
       a(n)=-11./23.
       val = 0.
       rhs(n)=24.*val/(23.*dx) +( -36./23.*f(n)+48*f(n-1)/23.-12.*f(n-2)/23.)/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
      end

      subroutine der2_s(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=0.5
       rhs(1)=(16.*val/5. -9./2*f(1)+f(2)+3./10*f(3))/(dx*dx)
       b(n)=1.
       c(n)=0.
       a(n)=0.5
       rhs(n)=(16.*val/5. -9./2*f(n)+f(n-1)+3./10*f(n-2))/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
c         do i=2,n-1
c         df(i)=(f(i+1)-2*f(i)+f(i-1))/(dx*dx)
c         enddo
c         df(1)=(-3*f(1)+f(2))/(dx*dx)
c         df(n)=(-3*f(n)+f(n-1))/(dx*dx)
      end

      subroutine der2_c(n,f,df,dx)
      integer n
      real f(n+1),df(n+1),a(n+1),b(n+1),c(n+1),d(n+1),rhs(n+1)

      do i=3,n-1
	 a(i)=2./11
	 b(i)=1.
         c(i)=2./11
         rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
      enddo
	b(2)=1.
        a(2)=1./10
        c(2)=1./10
        rhs(2) =6*(f(3)-2*f(2)+f(1))/(5*dx*dx)
	b(n)=1.
        a(n)=1./10
        c(n)=1./10
        rhs(n) =6*(f(n-1)-2*f(n)+f(n+1))/(5*dx*dx)
        a(1)=0
        b(1)=1
        c(1)=11. 
        rhs(1)=(13*f(1)-27*f(2)+15*f(3)-f(4))/(dx*dx)
        a(n+1)=11
        b(n+1)=1
        c(n+1)=0. 
        rhs(n+1)=(13*f(n+1)-27*f(n)+15*f(n-1)-f(n-2))/(dx*dx)
       call trid(n+1,a,b,c,rhs,d)
       do i=1,n+1	
         df(i)=rhs(i)
       enddo
      end



      subroutine interpol_c_s6(flag,N,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      RHS(i  )=2./3*(F(i)+F(i-1  ))


      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )


      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))*.5 + fac2*(F(i+1)+F(i-2))*.5
      enddo 
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      do i=3,n-2
	a1(i,i-1)=fac2*.5
	a1(i,i  )=fac1*.5
	a1(i,i+1)=fac1*.5
	a1(i,i+2)=fac2*.5
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      RHS(1) =(1./4)*f(0)+(3./2)*f(1)+f(2)/4
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      RHS(n) =(1./4)*f(n)+(3./2)*f(n-1)+f(n-2)/4
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
	  
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      
      subroutine interpol_c_s6_m(flag,N,NRHS,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N,NRHS)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      do j=1,nrhs
      RHS(i,j  )=2./3*(F(i,j)+F(i-1,j  ))
      enddo

      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )


      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)+F(i-1,j))*.5 + fac2*(F(i+1,j)+F(i-2,j))*.5
        enddo
      enddo 
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      do i=3,n-2
	a1(i,i-1)=fac2*.5
	a1(i,i  )=fac1*.5
	a1(i,i+1)=fac1*.5
	a1(i,i+2)=fac2*.5
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      do j=1,nrhs
      RHS(1,j) =(1./4)*f(0,j)+(3./2)*f(1,j)+f(2,j)/4
      enddo
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      do j=1,nrhs
      RHS(n,j) =(1./4)*f(n,j)+(3./2)*f(n-1,j)+f(n-2,j)/4
      enddo
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
	  
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call dgtsv(N,NRHS,A(2),B,C,RHS,N,info)
      do i=1,n
       do j=1,nrhs
       df(i,j)=rhs(i,j)
      enddo
      enddo
      end 
      
      subroutine der1_c_s6_m(flag,N,NRHS,F,dF,dx,t1,a1,r)
      real f(0:N,NRHS),df(N,NRHS),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N,NRHS)
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      do j=1,nrhs
      RHS(i ,j )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
      do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)-F(i-1,j))*dxi + fac2*(F(i+1,j)-F(i-2,j))*dxi /3.
      enddo
      enddo

      A(1)=0
      B(1)=1.
      C(1)=-1.
      do j=1,nrhs
      RHS(1,j) = dxi*(-f(0,j)+2*f(1,j)-f(2,j))
      enddo
      C(N)=0
      B(N)=1
      A(N)=-1.
      do j=1,nrhs
      RHS(N,j)=dxi*((f(n,j)-2*f(n-1,j)+f(n-2,j)))
      enddo
      call dgtsv(N,NRHS,A(2),b,c,rhs,n,info)
      do i=1,n
      do j=1,nrhs
       df(i,j)=rhs(i,j)
      enddo
      enddo
      end 
      


      subroutine der1_s_c6(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      do i=3,n-1
      a1(i,i-2)=-fac2*dxi/3
      a1(i,i-1)=-fac1*dxi
      a1(i,i  )= fac1*dxi
      a1(i,i+1)= fac2*dxi/3.
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  
      subroutine der1_s_c6_m(flag,N,nrhs,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,nrhs),dF(0:N,nrhs),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,nrhs),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i ,j )=3*(3-1./11)/8.*(F(i,j)-F(i-1 ,j ))*dxi
      enddo
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)-F(i-1,j))*dxi + fac2*(F(i+1,j)-F(i-2,j))*dxi /3.
        enddo
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      do i=3,n-1
      a1(i,i-2)=-fac2*dxi/3
      a1(i,i-1)=-fac1*dxi
      a1(i,i  )= fac1*dxi
      a1(i,i+1)= fac2*dxi/3.
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3,j  ))*dxi
      enddo
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.   
      do j=1,nrhs
      RHS(N+1 ,j )=(-ac*F(N,j)-bc*F(N-1,j)-cc*F(N-2,j))*dxi
      enddo
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,N+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end
  
      subroutine derw_s_c(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo
     

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  
      subroutine interpol_s_c6(flag,N,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k
      character*1 flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      t1 =0.
      a1 =0
      r  =0
      do i=2,N
	A   (i) = 1./6
        B   (i) = 1.
        C   (i) = 1./6
        RHS (i) = 2./3*(F(i)+F(i-1))
      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))*.5 + fac2*(F(i+1)+F(i-2))*.5
      enddo

	 
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=2./3
      a1(i,i  )=2./3
      enddo
      do i=3,n-1
      a1(i,i-2)=fac2*.5
      a1(i,i-1)=fac1*.5
      a1(i,i  )=fac1*.5
      a1(i,i+1)=fac2*.5
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(ac*F(N)+bc*F(N-1)+cc*F(N-2))
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end

      subroutine interpol_s_c6_m(flag,N,nrhs,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k,info,nrhs
      character*1 flag
      real F(N,nrhs),dF(0:N,nrhs),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,nrhs),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      t1 =0.
      a1 =0
      r  =0
      write(*,*)'dentro interpol_s_c6_m'
      do i=2,N
	A   (i) = 1./6
        B   (i) = 1.
        C   (i) = 1./6
        do j=1,nrhs
        RHS (i,j) = 2./3*(F(i,j)+F(i-1,j))
        enddo
      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)+F(i-1,j))*.5 + fac2*(F(i+1,j)+F(i-2,j))*.5
        enddo
      enddo

	 
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=2./3
      a1(i,i  )=2./3
      enddo
      do i=3,n-1
      a1(i,i-2)=fac2*.5
      a1(i,i-1)=fac1*.5
      a1(i,i  )=fac1*.5
      a1(i,i+1)=fac2*.5
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      do j=1,nrhs
      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3 ,j ))
      enddo
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
      RHS(N+1 ,j )=(ac*F(N,j)+bc*F(N-1,j)+cc*F(N-2,j))
      enddo
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
       enddo
      enddo
      end

      subroutine interpol_c_s4(flag,N,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      RHS(i  )=2./3*(F(i)+F(i-1  ))
      enddo
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      RHS(1) =(1./4)*f(0)+(3./2)*f(1)+f(2)/4
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      RHS(n) =(1./4)*f(n)+(3./2)*f(n-1)+f(n-2)/4
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      
      subroutine der1_c_s4(flag,N,F,dF,dx,t1,a1,r)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=-3*(3-1./11)/8.*dxi
      a1(i,i+1)=3*(3-1./11)/8.*dxi
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=-1.
      RHS(1) = dxi*(-f(0)+2*f(1)-f(2))
      C(N)=0
      B(N)=1
      A(N)=-1.
      RHS(N)=dxi*((f(n)-2*f(n-1)+f(n-2)))
      if (flag.eq.'m') then
      a1(1,1)=-dxi
      a1(1,2)=2*dxi
      a1(1,3)=-dxi
      a1(n,n+1)=dxi
      a1(n,n  )=-2*dxi
      a1(n,n-1)=dxi
       do i=2,n-1
	  t1(i,i)  =b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td=0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      


      subroutine der1_s_c4(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  

      subroutine interpol_s_c4(flag,N,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k
      character*1 flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      alpha =1./6
      fac1=2./3.
      t1 =0.
      a1 =0
      r  =0
      do i=2,N
	A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))
      enddo
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=fac1
      a1(i,i  )=fac1
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(ac*F(N)+bc*F(N-1)+cc*F(N-2))
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=900)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END



      subroutine der1w_s_c6(N,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi

      A  (N+1  )=15
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(15*f(n)-50*f(n-1)/3+0.6*f(n-2))*dxi
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end
  
      subroutine der1w_s_c6_m(N,NRHS,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
c      RHS(1 ,j )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3,j  ))*dxi
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j  )/3-0.6*F(3,j ))*dxi
      enddo
      A  (N+1  )=15
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
      RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j))*dxi
      enddo
c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end
 

      subroutine der1c_s_c6_m(N,NRHS,F,dF,dx)
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real Qval,mr_c(0:N+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15!23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
c      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j)+cc*F(3,j))*dxi
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j)-16./15.*1.0 )*dxi
      enddo
      A  (N+1  )=15.
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
       RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j)+16./15.*1.0)*dxi
      enddo

c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo

      end


 
!     subroutine der1c_s_c6_m(N,NRHS,F,dF,dx,Twh,Twc)
!     implicit none
!     integer i,n,j,k,info,nrhs
!     character flag
!     real F(N,NRHS),dF(0:N,NRHS),dx
!     real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
!     real dxi,ac,bc,cc,alpha,fac1,fac2
!     real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
!     real Qval,mr_c(0:N+1),Twh,Twc
!     a1=0
!     t1=0
!     r=0
!     dxi = 1./dx
!     alpha= 9./62
!     fac1 = (3./8)*(3-2*alpha)
!     fac2 = (1./8)*(-1+22*alpha)
!     do i=2,n
!     A  (i  )=1./22
!     B  (i  )=1
!     C  (i  )=1./22
!     do j=1,nrhs
!     RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
!     enddo
!     enddo
!     A  (1  )=0.
!     B  (1  )=1.
!     C  (1  )=15!23.
!     ac =-(23**2+71)/24.
!     bc = (7*23+47)/8.
!     cc = (23-31)/8.
!     do j=1,nrhs
!     RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j)-16./15.*Twh)*dxi
!     enddo
!     A  (N+1  )=15.
!     B  (N+1  )=1.
!     C  (N+1  )=0.
!     do j=1,nrhs
!      RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j)+16./15.*Twc)*dxi
!     enddo
!

!     call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
!     do i=1,N+1
!      do j=1,nrhs
!      dF(i-1,j)=RHS(i,j)
!     enddo
!     enddo
!     end
! 
     
      subroutine der1_c_s6(flag,N,F,dF,dx,t1,a1,r)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N)
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      dxi=1./dx
      i=2
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      i=n-1
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo

      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=-3*(3-1./11)/8.*dxi
      a1(i,i+1)=3*(3-1./11)/8.*dxi
      enddo

      do i=3,n-2
      a1(i,i-1)=-fac2*dxi/3.
      a1(i,i  )=-fac1*dxi
      a1(i,i+1)=fac1*dxi
      a1(i,i+2)=fac2*dxi/3.
      enddo
      endif
      A(1)=0
      B(1)=1.
      C(1)=-1.
      RHS(1) = dxi*(-f(0)+2*f(1)-f(2))
      C(N)=0
      B(N)=1
      A(N)=-1.
      RHS(N)=dxi*((f(n)-2*f(n-1)+f(n-2)))
      if (flag.eq.'m') then
      a1(1,1)=-dxi
      a1(1,2)=2*dxi
      a1(1,3)=-dxi
      a1(n,n+1)=dxi
      a1(n,n  )=-2*dxi
      a1(n,n-1)=dxi
       do i=2,n-1
          t1(i,i)  =b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td=0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
         do k=1,n+1
         do i=1,n
          r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end

     

 
      function f_colebrook(Re,bulk)
      real cole
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - cole(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_colebrook=f
      end
      function cole(f,Reb)
       cole = ( 1/sqrt(f) +2.*log10(2.51/(Reb*sqrt(f))))
      end

      function f_shokling(Re,bulk)
      real shok 
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - shok(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_shokling=f
      end
      function shok(f,Reb)
       shok= ( 1/sqrt(f) -1.93*log10(Reb*sqrt(f))+0.537)
      end

      subroutine output_23D(istap,rank)
      include 'par_post.txt'
      include 'common.txt'
      integer istap
      character*5 cha
      character*5 cha2
      real um(imax),vm(imax),wm(imax)
      real umm(imax),vmm(imax),wmm(imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),hi_s(imax),hi_c(0:imax)

      real fl(0:i1,jmax/p_row,kmax/p_col) , xx(3),yy(3)
      real vor  (0:i1,jmax/p_row,kmax/p_col)
      real vot  (0:i1,jmax/p_row,kmax/p_col)
      real voz  (0:i1,jmax/p_row,kmax/p_col)
      real votot(0:i1,jmax/p_row,kmax/p_col)
     
      
      integer idex,jdex,kdex,i_dex,j_dex,k_dex


      kdex  = ( nrank - (nrank/p_col)*p_col) * (kmax/p_col)
      jdex  = ((nrank/p_col)) * (jmax/p_row)
 

      call cnvstr(rank,cha )
      call cnvstr(rank,cha2)
     
      call cnvstr(istap,cha)
      call cnvstr(rank,cha2)



      open(45,file ='tec_rt3d.'//cha2)
      write(45,*) ' VARIABLES ="X", "Y","Z", "U-vel","V-vel","W-vel","P", "L2", "VR" ,"VT", "VZ" ,"VTOT" '
      write(45,*) ' ZONE I= ',imax, ', J= ', jmax/p_row -2 , 'k = ',kmax/p_col -2 ,'  F=POINT '
      do k=2,kmax/p_col-1
      do j=2,jmax/p_row-1
	do i=1,imax
	  write(45,'(14E16.7)') rp(i)*cos((j+jdex)*dtheta),rp(i)*sin((j+jdex)*dtheta),(k+kdex)*dz,
     ^ unew(i,j,k),vnew(i,j,k),wnew(i,j,k),p(i,j,k),fl(i,j,k),vor(i,j,k),vot(i,j,k),voz(i,j,k),votot(i,j,k)
	enddo
      enddo
      enddo
      close(45)
      end



C
C      ________________________________________________________
C     |                                                        |
C     |            SORT AN ARRAY IN INCREASING ORDER           |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         X     --ARRAY OF NUMBERS                       |
C     |                                                        |
C     |         Y     --WORKING ARRAY (LENGTH  AT LEAST N)     |
C     |                                                        |
C     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         X     --SORTED ARRAY                           |
C     |________________________________________________________|
C
      SUBROUTINE SORT(X,Y,N)
      REAL X(1),Y(1),S,T
      INTEGER I,J,K,L,M,N
      I = 1
10    K = I
20    J = I
      I = I + 1
      IF ( J .EQ. N ) GOTO 30
      IF ( X(I) .GE. X(J) ) GOTO 20
      Y(K) = I
      GOTO 10
30    IF ( K .EQ. 1 ) RETURN
      Y(K) = N + 1
40    M = 1
      L = 1
50    I = L
      IF ( I .GT. N ) GOTO 120
      S = X(I)
      J = Y(I)
      K = J
      IF ( J .GT. N ) GOTO 100
      T = X(J)
      L = Y(J)
      X(I) = L
60    IF ( S .GT. T ) GOTO 70
      Y(M) = S
      M = M + 1
      I = I + 1
      IF ( I .EQ. K ) GOTO 80
      S = X(I)
      GOTO 60
70    Y(M)= T
      M = M + 1
      J = J + 1
      IF ( J .EQ. L ) GOTO 110
      T = X(J)
      GOTO 60
80    Y(M) = T
      K = M + L - J
      I = J - M
90    M = M + 1
      IF ( M .EQ. K ) GOTO 50
      Y(M) = X(M+I)
      GOTO 90
100   X(I) = J
      L = J
110   Y(M) = S
      K = M + K - I
      I = I - M
      GOTO 90
120   I = 1
130   K = I
      J = X(I)
140   X(I) = Y(I)
      I = I + 1
      IF ( I .LT. J ) GOTO 140
      Y(K) = I
      IF ( I .LE. N ) GOTO 130
      IF ( K .EQ. 1 ) RETURN
      GOTO 40
      END
