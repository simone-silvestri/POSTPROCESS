c234567
      program mpialltoall_test
      include 'mpif.h'
      parameter (nx=20,ny=20,nz=20,px=5)
      real*8 Ax (nx,ny/px,nz/px)
      real*8 Ay (nx/px,ny,nz/px)
      real*8 Az (nx/px,ny/Px,nz)

  
      SUBROUTINE init_transpose
        INCLUDE 'param.txt'
C ...  Globals
	integer Xii(Nx),Xkk(Nx,Mt)
	common /XPOSE/ Xii,Xkk
	integer Tkk(Nt),Tii(Nt,Mx)
	common /TPOSE/ Tkk,Tii

C ...  Locals
	  
	do k = 1,Mt
	  do i = 1,Nx
	    Xii(i) = MOD(i-1,Mx)+1
	    Xkk(i,k) = INT((i-1)/Mx)*Mt+k
	  end do
	end do
	
	do i = 1,Mx
	  do k = 1,Nt
	    Tkk(k) = MOD(k-1,Mt) + 1
	    Tii(k,i) = INT((k-1)/Mt)*Mx + i
	  end do
	end do
	  
      RETURN
      END


C ... Perfroms data transform across the nodes.  The data
C     starts out Q(j,k,i) node distributed in x and is transformed
C     to Qt(i,j,k) node distributed in theta.

      SUBROUTINE t_j2k(Q,Qt,rank)
        INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
	real  Q(Nr,Nt,Mx),Qt(Nr,Mt,Nx)
		
C ...  Locals
        
C ...  Globals
	integer Xii(Nx),Xkk(Nx,Mt)
	common /XPOSE/ Xii,Xkk
	
	real*8  W1(Mx,Nr,Nt),W2(Mx,Nr,Nt)
	  
	  do j = 1,Nr
	  do k = 1,Nt
	  do i = 1,Mx
		W1(i,j,k) = Q(j,k,i)
	      end do
	    end do
	  end do
          
	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
	  
	  do k = 1,Mt
	    do j = 1,Nr
	      do i = 1,Nx
		Qt(j,k,i) = W2(Xii(i),j,Xkk(i,k))
	      end do
	    end do
	  end do
	  
		
      RETURN
      END
      
      
C ... Performs data transform across the nodes.  The data
C     starts out (i,j,k) node distributed in theta and is transformed
C     to (j,k,i) node distributed in x.

      SUBROUTINE t_k2j(Q,Qt,rank)
        INCLUDE 'param.txt'
        INCLUDE 'mpif.h'

	real*8 Q(Nr,Mt,Nx),Qt(Nr,Nt,Mx)

C ...  Locals

C ...  Globals
	integer Tkk(Nt),Tii(Nt,Mx)
	common /TPOSE/ Tkk,Tii
	  	  
	real W2(Nr,Mt,Nx)
	  
	  call MPI_ALLTOALL(Q,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
     
	  do i = 1,Mx
	    do k = 1,Nt
	      do j = 1,Nr
		Qt(j,k,i) = W2(j,Tkk(k),Tii(k,i))
	      end do
	    end do
	  end do
	  

      RETURN
      END

