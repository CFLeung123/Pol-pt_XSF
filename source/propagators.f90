


	MODULE Parameters
	!
	! This module contains:
	!
	! 		- Definitions of essential parameters for the rest of the calculation.
	!		- Definition of several utility functions and subroutines
	!		- Definition of traces of products of 4x4 matrices
	!		- Definition of functions for naming correlation functions
	!
	implicit none
	integer :: c0, c1	! Parameters such that x_0 in [c0,c1]
	complex*16,parameter :: imag=(0.0d0,1.0d0)
	real(kind=8),parameter :: pi=3.1415926535897932384626433832d0
	real(kind=8),parameter :: dimrep = 3.0d0
	real(kind=8),parameter :: C2_R = 4.0d0/3.0d0
	integer, dimension(2,2) :: f1f2 !flavour keep_tracker
	integer,dimension(:,:,:),allocatable :: multiplicity

	CONTAINS

	SUBROUTINE flav_counter_init()
	!
	! Sets a value to index flavour pairs.
	!
	f1f2(1,1)=1	!uu
	f1f2(1,2)=2	!ud
	f1f2(2,1)=3	!du
	f1f2(2,2)=4	!dd

	END SUBROUTINE flav_counter_init




	SUBROUTINE momentum_degeneracy(l,mom_deg)
	!
	! Computes the degeneracy of a momentum vector.
	!
	integer,intent(in) :: l, mom_deg
	integer :: i,j,k
	integer :: nequals

	allocate(multiplicity(0:l-1,0:l-1,0:l-1))

	multiplicity(:,:,:) = 0
	nequals=0


	if(mom_deg.eq.1)then

	do i=0,l-1
	do j=i,l-1
	do k=j,l-1

		nequals=0
		if(i.eq.j)then
	  	nequals=2
		  if(j.eq.k)then
		    nequals=3
		  else
		  endif
		else
		  if(i.eq.k)then
		    nequals=2
		  else
		    if(j.eq.k)then
		      nequals=2
		    else
		    endif
		  endif
		endif

		if(nequals.eq.0)then
		  multiplicity(i,j,k)=6
		else
		  if(nequals.eq.2)then
		    multiplicity(i,j,k)=3
		  else
		    multiplicity(i,j,k)=1
		  endif
		endif

	enddo
	enddo
	enddo

	else

		multiplicity(:,:,:)=1

	endif

	END SUBROUTINE momentum_degeneracy





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	SUBROUTINE adjoint_matrix(A,B)
	!
	! Subroutine to calculate the adjoint of a 4x4 matrix
	!
	implicit none

	complex(kind=8),dimension(4,4),intent(in) :: A
	complex(kind=8),dimension(4,4),intent(out) :: B
	integer :: i,j

	do i=1,4
	do j=1,4

	B(j,i) = dcmplx(real(A(i,j)),-aimag(A(i,j)))

	enddo
	enddo

	END SUBROUTINE adjoint_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE TemporalExtent(l,deriv)
	!
	! Sets the limits of the temporal extent in which the different correlation
	! functions are evaluated.
	!
	implicit none
	integer :: l
	integer :: deriv

	if(even_or_odd(l).eq.1)then
	  c0 = l/2-1
	  c1 = l/2+1
	  if(deriv .eq. 0)then		!If we don't want derivatives of the correlation functions
	   c0 = c0 + 1				!this will fix the calculation to be done only at x0=l/2.
	   c1 = c0
	  endif
	else
	  c0=(l-1)/2-1
	  c1=(l+1)/2+1
	  if(deriv .eq. 0)then		!If we don't want derivatives of the correlation functions
	   c0 = c0 + 1				!this will fix the calculation to be done only at x0=(l-1)/2
	   c1 = c1-1				!and at x0=(l+1)/2.
	  endif
	endif


	END SUBROUTINE TemporalExtent



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE loop_limits(l,s0,t0_a,t0_b,u0_a,u0_b)
!	-This subroutine defines correctly, for the XSF
!	 the temporal indeces inside temporal
!	 sums when evaluating Sigma1, Sigma2, Gamma1 and the loop diagrams.
!	-It takes avantage of the delta functions in the vertex functions.
!	-NOTE that it is adapted to work only with the intel compiler,
!	 since the boolean values are:
!			true=-1,    false=0
!	-For a given temporal index s0, it returns the range over which the other
!	 indeces t0 and u0 should be integrated, which are always adjacent points.
!	-If s0 is close to the boundaries (1 or l-1), then the limits of the ranges
!	 are modified to avoid summing outside the SF.
!
!	-Input: -s0
!
!	-Output: -t0_a, t0_b: range for t0
!		 -x0_a, x0_b: range for x0
	implicit none

	integer,intent(in) :: l,s0
	integer,intent(out) :: t0_a,t0_b,u0_a,u0_b
	integer :: c,d,e

		c = merge(-1, 0, s0 .eq. 0)
		d = merge(-1, 0, s0 .eq. l)
		e = merge(-1, 0, s0 .eq. (l-1))

!	   if(s0.eq.1)then
!	    c=-1
!	   else
!	    c=0
!	   endif
!	   if(s0.eq.(l-1))then
!	    d=-1
!	   else
!	    d=0
!	   endif

!	    t0_a = s0 - 1 - c
!	    t0_b = s0 + 1 + d
!	    u0_a = s0 - 1
!	    u0_b = s0 + 1 + d

	    t0_a = s0 - 1 - c
	    t0_b = s0 + 1 + d
	    u0_a = s0 - 1 - c
	    u0_b = s0 + 1 + e + 2*d


	END SUBROUTINE loop_limits


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER FUNCTION  n_limit(deg,n)
	!
	! Function to define the lower limits of a momentum loop deppending on wether
	! we are taking advantage of momentum degeneracy or not.
	!
	integer :: deg, n

	if(deg.eq.1)then
	  n_limit = n
	else
	  n_limit = 0
	endif

	END FUNCTION n_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	INTEGER FUNCTION m_period(n,l)
	!
	! Function that implements periodic boundary conditions.
	!
	integer :: n,l

	if (n == 0) then
	  m_period = 0
	else
	  m_period = l - n
	endif

	END FUNCTION m_period



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	REAL(kind=8) FUNCTION delta(a,b)
	!
	! Kronecker delta function
	!
	integer:: a,b

	if(a.eq.b)then
	  delta=1.0d0
	else
	  delta=0.0d0
	endif

	END FUNCTION delta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	REAL(kind=8) FUNCTION omega(a,b)
	!
	! Heavyside function
	!
	integer :: a,b

	if(a.gt.b)then
	  omega=1.0d0
	else
	  omega=0.0d0
	endif

	END FUNCTION omega

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	INTEGER Function even_or_odd(l)
	!
	! Function to determine wether the number l is even or odd
	!
	integer :: l
	real :: a

	a = (-1.0)**l
	if(a.gt.(0.0))then
	  even_or_odd = 1	!l is even number
	else
	  even_or_odd = 0	!l is odd number
	endif

	END FUNCTION even_or_odd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	COMPLEX(kind=8) FUNCTION Trace2(A,B)
	!
	! Function to calculate the trace of a product of 2 4x4 matrices:
	!
	complex(kind=8),dimension(4,4) :: A,B
	integer :: i,j

	Trace2 = dcmplx(0.0d0,0.0d0)

	do i=1,4
	do j=1,4

	  Trace2 = Trace2 + A(j,i)*B(i,j)

	enddo
	enddo

	END FUNCTION Trace2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	COMPLEX(kind=8) FUNCTION Trace3(A,B,C)
	!
	! Function to calculate the trace of a product of 3 4x4 matrices:
	!
	complex(kind=8),dimension(4,4) :: A,B,C
	complex(kind=8) :: Acum1
	integer :: i,j,k

	Trace3 = dcmplx(0.0d0,0.0d0)

	do i=1,4
	do k=1,4

	  Acum1 = dcmplx(0.0d0,0.0d0)

	  do j=1,4

 	    Acum1 = Acum1 + A(i,j)*B(j,k)

	  enddo

          Trace3 = Trace3 + Acum1*C(k,i)

	enddo
	enddo

	END FUNCTION Trace3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	COMPLEX(kind=8) FUNCTION Trace4(A,B,C,D)
	!
	! Function to calculate the trace of a product of 4 4x4 matrices:
	!

	complex(kind=8),dimension(4,4) :: A,B,C,D
	complex(kind=8) :: Acum1, Acum2
	integer :: i,j,k,l

	Trace4 = dcmplx(0.0d0,0.0d0)

	do i=1,4
	do l=1,4

	  Acum1 = dcmplx(0.0d0,0.0d0)

  	  do k=1,4

	    Acum2 = dcmplx(0.0d0,0.0d0)

 	    do j=1,4

	      Acum2 = Acum2 + A(i,j)*B(j,k)

	    enddo

	    Acum1 = Acum1 + Acum2*C(k,l)

	  enddo

	  Trace4 = Trace4 + Acum1*D(l,i)

	enddo
	enddo

	END FUNCTION Trace4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	COMPLEX(kind=8) FUNCTION Trace5(A,B,C,D,E)
	!
	! Function to calculate the trace of a product of 5 4x4 matrices:
	!

	complex(kind=8),dimension(4,4) :: A,B,C,D,E
	complex(kind=8) :: Acum1, Acum2, Acum3
	integer :: i,j,k,l,m

	Trace5 = dcmplx(0.0d0,0.0d0)

	do i=1,4
	do m=1,4

	  Acum1 = dcmplx(0.0d0,0.0d0)

  	  do l=1,4

	    Acum2 = dcmplx(0.0d0,0.0d0)

 	    do k=1,4

	      Acum3 = dcmplx(0.0d0,0.0d0)

	        do j=1,4

		  Acum3 = Acum3 + A(i,j)*B(j,k)

		enddo

		Acum2 = Acum2 + Acum3*C(k,l)

	    enddo

	    Acum1 = Acum1 + Acum2*D(l,m)

	  enddo

	  Trace5 = Trace5 + Acum1*E(m,i)

	enddo
	enddo

	END FUNCTION Trace5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	COMPLEX(kind=8) FUNCTION Trace6(A,B,C,D,E,F)
	!
	! Function to calculate the trace of a product of 6 4x4 matrices:
	!
	complex(kind=8),dimension(4,4) :: A,B,C,D,E,F
	complex(kind=8) :: Acum1, Acum2, Acum3, Acum4
	integer :: i,j,k,l,m,n


	Trace6 = dcmplx(0.0d0,0.0d0)

	do i=1,4
	do n=1,4

	  Acum1 = dcmplx(0.0d0,0.0d0)

  	  do m=1,4

	    Acum2 = dcmplx(0.0d0,0.0d0)

 	    do l=1,4

	      Acum3 = dcmplx(0.0d0,0.0d0)

	      do k=1,4

		Acum4 = dcmplx(0.0d0,0.0d0)

		do j=1,4

		  Acum4 = Acum4 + A(i,j)*B(j,k)

		enddo

		Acum3 = Acum3 + Acum4*C(k,l)

	      enddo

	      Acum2 = Acum2 + Acum3*D(l,m)

	    enddo

	    Acum1 = Acum1 + Acum2*E(m,n)

	  enddo

	  Trace6 = Trace6 + Acum1*F(n,i)

	enddo
	enddo

	END FUNCTION Trace6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	CHARACTER*2 FUNCTION  corrname_gXff(i)
	!
	! Returns a string with the name of the gXff correlation functions
	!
	integer :: i

	if(i.eq.1)then
	  corrname_gXff = 'gA'
	else
	endif

	if(i.eq.2)then
	  corrname_gXff = 'gP'
	else
	endif

	if(i.eq.3)then
	  corrname_gXff = 'gV'
	else
	endif

	if(i.eq.4)then
	  corrname_gXff = 'gS'
	else
	endif

	END FUNCTION corrname_gXff




	CHARACTER*3 FUNCTION  corrname_lYff(i)
	!
	! Returns a string with the name of the lYff correlation functions
	!
	integer :: i

	if(i.eq.1)then
	  corrname_lYff = 'lA_'
	else
	endif

	if(i.eq.2)then
	  corrname_lYff = 'lV_'
	else
	endif

	if(i.eq.3)then
	  corrname_lYff = 'lT_'
	else
	endif

	if(i.eq.4)then
	  corrname_lYff = 'lTt'
	else
	endif

	END FUNCTION corrname_lYff



	CHARACTER*2 FUNCTION  flavname(i)
	!
	! Returns a string with the flavour combinations uu, ud, du, dd.
	!
	integer :: i

	if(i.eq.1)then
	  flavname = 'uu'
	else
	endif

	if(i.eq.2)then
	  flavname = 'ud'
	else
	endif

	if(i.eq.3)then
	  flavname = 'du'
	else
	endif

	if(i.eq.4)then
	  flavname = 'dd'
	else
	endif

	END FUNCTION flavname


	character(len=4) function str2(k)
	!
	! Convert an integer to string.
	!
    	integer, intent(in) :: k
    	write (str2, *) k
    	str2 = adjustl(str2)
	end function str2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	CHARACTER*5 FUNCTION  corrname_Bilff(i)
	!
	! Returns a string with the name of the gXff correlation functions
	!
	integer :: i

	if(i.eq.1)then
	  corrname_Bilff = 'gA0'
	else
	endif

	if(i.eq.2)then
	  corrname_Bilff = 'gP'
	else
	endif

	if(i.eq.3)then
	  corrname_Bilff = 'gV0'
	else
	endif

	if(i.eq.4)then
	  corrname_Bilff = 'gS'
	else
	endif


	if(i.eq.5)then
	  corrname_Bilff = 'gA1'
	else
	endif

	if(i.eq.6)then
	  corrname_Bilff = 'gA2'
	else
	endif

	if(i.eq.7)then
	  corrname_Bilff = 'gA3'
	else
	endif

	if(i.eq.8)then
	  corrname_Bilff = 'gV1'
	else
	endif

	if(i.eq.9)then
	  corrname_Bilff = 'gV2'
	else
	endif

	if(i.eq.10)then
	  corrname_Bilff = 'gV3'
	else
	endif

	if(i.eq.11)then
	  corrname_Bilff = 'gT10'
	else
	endif

	if(i.eq.12)then
	  corrname_Bilff = 'gT20'
	else
	endif

	if(i.eq.13)then
	  corrname_Bilff = 'gT30'
	else
	endif

	if(i.eq.14)then
	  corrname_Bilff = 'gTt10'
	else
	endif

	if(i.eq.15)then
	  corrname_Bilff = 'gTt20'
	else
	endif

	if(i.eq.16)then
	  corrname_Bilff = 'gTt30'
	else
	endif

	END FUNCTION corrname_Bilff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	CHARACTER*6 FUNCTION  corrname_4f(i)
	!
	! Returns a string with the name of a correlation functions of a 4 fermion
	! operator.
	!
	integer :: i

	if(i.eq.1)then
	  corrname_4f = 'II____'
	else
	endif

	if(i.eq.2)then
	  corrname_4f = 'g5____'
	else
	endif

	if(i.eq.3)then
	  corrname_4f = 'g0____'
	else
	endif

	if(i.eq.4)then
	  corrname_4f = 'g1____'
	else
	endif


	if(i.eq.5)then
	  corrname_4f = 'g2____'
	else
	endif

	if(i.eq.6)then
	  corrname_4f = 'g3____'
	else
	endif

	if(i.eq.7)then
	  corrname_4f = 'g0g5__'
	else
	endif

	if(i.eq.8)then
	  corrname_4f = 'g1g5__'
	else
	endif

	if(i.eq.9)then
	  corrname_4f = 'g2g5__'
	else
	endif

	if(i.eq.10)then
	  corrname_4f = 'g3g5__'
	else
	endif

	if(i.eq.11)then
	  corrname_4f = 'is00__'
	else
	endif

	if(i.eq.12)then
	  corrname_4f = 'is01__'
	else
	endif

	if(i.eq.13)then
	  corrname_4f = 'is02__'
	else
	endif

	if(i.eq.14)then
	  corrname_4f = 'is03__'
	else
	endif

	if(i.eq.15)then
	  corrname_4f = 'is10__'
	else
	endif

	if(i.eq.16)then
	  corrname_4f = 'is11__'
	else
	endif

	if(i.eq.17)then
	  corrname_4f = 'is12__'
	else
	endif

	if(i.eq.18)then
	  corrname_4f = 'is13__'
	else
	endif

	if(i.eq.19)then
	  corrname_4f = 'is20__'
	else
	endif

	if(i.eq.20)then
	  corrname_4f = 'is21__'
	else
	endif

	if(i.eq.21)then
	  corrname_4f = 'is22__'
	else
	endif

	if(i.eq.22)then
	  corrname_4f = 'is23__'
	else
	endif

	if(i.eq.23)then
	  corrname_4f = 'is30__'
	else
	endif

	if(i.eq.24)then
	  corrname_4f = 'is31__'
	else
	endif

	if(i.eq.25)then
	  corrname_4f = 'is32__'
	else
	endif

	if(i.eq.26)then
	  corrname_4f = 'is33__'
	else
	endif

	if(i.eq.27)then
	  corrname_4f = 'ig5s00'
	else
	endif

	if(i.eq.28)then
	  corrname_4f = 'ig5s01'
	else
	endif

	if(i.eq.29)then
	  corrname_4f = 'ig5s02'
	else
	endif

	if(i.eq.30)then
	  corrname_4f = 'ig5s03'
	else
	endif

	if(i.eq.31)then
	  corrname_4f = 'ig5s10'
	else
	endif

	if(i.eq.32)then
	  corrname_4f = 'ig5s11'
	else
	endif

	if(i.eq.33)then
	  corrname_4f = 'ig5s12'
	else
	endif

	if(i.eq.34)then
	  corrname_4f = 'ig5s13'
	else
	endif

	if(i.eq.35)then
	  corrname_4f = 'ig5s20'
	else
	endif

	if(i.eq.36)then
	  corrname_4f = 'ig5s21'
	else
	endif

	if(i.eq.37)then
	  corrname_4f = 'ig5s22'
	else
	endif

	if(i.eq.38)then
	  corrname_4f = 'ig5s23'
	else
	endif

	if(i.eq.39)then
	  corrname_4f = 'ig5s30'
	else
	endif

	if(i.eq.40)then
	  corrname_4f = 'ig5s31'
	else
	endif

	if(i.eq.41)then
	  corrname_4f = 'ig5s32'
	else
	endif

	if(i.eq.42)then
	  corrname_4f = 'ig5s33'
	else
	endif


	END FUNCTION corrname_4f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	END MODULE Parameters
