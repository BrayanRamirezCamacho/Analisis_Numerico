!---------------------------------------------------------
!Este programa calcula el epsilon de una computadora.
!
!---------------------------------------------------------

PROGRAM epsilon

IMPLICIT NONE
REAL(8)::e
real(16)::epsi
INTEGER::i,n
print*,"Quiere calcular epsilon con presicion sencilla o doble?" 
print*, "Sencilla: 1"
print*, "Doble: 2"

read*,n

IF (n==1) then	
	e=1
	DO 
	e=e/(2.)
	IF (1+e==1) EXIT
	END DO
	PRINT*, "epsilon=",e

else if (n==2)then
	epsi=1
	DO 
	epsi=epsi/(2.)
	IF (1+epsi==1) EXIT
	END DO
	PRINT*, "epsilon=",epsi

else

	PRINT*,"Introduzca una opcion valida."


end if

END PROGRAM 
