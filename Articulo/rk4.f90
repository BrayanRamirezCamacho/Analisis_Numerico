  !!El objetivo de este programa es el de aproximar la solución de una
  !!ecuación diferencial (E.D) utilizado el método de Runge-Kutta de
  !!orden 4 (RK4).
  !!-------------------------------------------------------------------

PROGRAM rungekuta4

IMPLICIT NONE
!Diccionario de variables
 
REAL,PARAMETER::M=3,K=50,L=14
REAL::V0,X01,X02,h,TF,T0
INTEGER::n

PRINT*,"introduzca ambas posiciones iniciales"
READ*,X01,X02
PRINT*,"introduzca el tiempo"
READ*,TF
print*,"Introduzca n"
read*, n

T0=0
!n=50
V0=0
h=(TF-T0)/n

call SRK4(T0,X01,X02,h,n)

end program rungekuta4

!_______________________________________________________________
SUBROUTINE SRK4(T0,X01,X02,h,n)
IMPLICIT NONE
                                        
REAL,INTENT(IN)::T0,X01,X02,h
REAL::V0
INTEGER,INTENT(IN)::n
REAL::X1,T,V1,V2,X2,XA1,XA2
REAL,DIMENSION(1:4,1:4)::K
INTEGER::i
REAL,EXTERNAL::ED1

X1=X01
X2=X02
V1=V0
V2=V0
T=T0
XA1=X1
XA2=X2

OPEN(UNIT=1,FILE="Datos",STATUS="UNKNOWN",ACTION="WRITE")

WRITE(1,*)T,",",X1,",",X2,",",V1,",",V2

DO i=1,n
	K(1,1)=V1
	K(2,1)=ED1(X1,XA1)
	K(3,1)=V2
	K(4,1)=ED1(X2,XA2)
	
	K(1,2)=V1+ h*K(2,1)/2.
	K(2,2)=ED1(X1+h*K(1,1)/2.,XA1)
	K(3,2)=V2+h*K(4,1)/2.
	K(4,2)=ED1(X2+ h*K(3,1)/2.,XA2)
	
	K(1,3)=V1+ h*K(2,2)/2.
	K(2,3)=ED1(X1+h*K(1,2)/2.,XA1)
	K(3,3)=V2+h*K(4,2)/2.
	K(4,3)=ED1(X2+ h*K(3,2)/2.,XA2)
	
	K(1,4)=V1+ h*K(2,3)/2.
	K(2,4)=ED1(X1+h*K(1,3)/2.,XA1)
	K(3,4)=V2+h*K(4,3)/2.
	K(4,4)=ED1(X2+ h*K(3,3)/2.,XA2)

	X1=X1 + h/6.*(K(1,1) + 2*(K(1,2) + K(1,3)) + K(1,4))
	V1=V1 + h/6.*(K(2,1) + 2*(K(2,2) + K(2,3)) + K(2,4))
	X2=X2 + h/6.*(K(3,1) + 2*(K(3,2) + K(3,3)) + K(3,4))
	V2=V2 + h/6.*(K(4,1) + 2*(K(4,2) + K(4,3)) + K(4,4))	
	XA1=X1
	XA2=X2
	T=T + h


WRITE(1,*)T,",",X1,",",X2,",",V1,",",V2

END DO

 CLOSE(1)

END SUBROUTINE SRK4

!________________________________________________


FUNCTION ED1(X,Y)
IMPLICIT NONE
REAL::ED1,X,Y

ED1=(-1)*(100/3.)*X + (50/3.)*Y

END FUNCTION ED1



























