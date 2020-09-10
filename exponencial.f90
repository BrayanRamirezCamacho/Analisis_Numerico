!------------------------------------------------------------------
!Este programa evalua la funcion f(x)=exp(x)-exp(-2x), y la 
!funcion f(x)=x-sen(x)
!usando la serie de MacLaurin, tomando n terminos de la sumatoria
! y evaluando la funcion en un valor de x; ambos escogidos por el
!usuario.
!------------------------------------------------------------------


PROGRAM exponencial

IMPLICIT NONE
	REAL::x,f_1,f_2
	INTEGER::i,n
	INTEGER,EXTERNAL::factorial
	INTEGER::p

PRINT*,"Que funcion desea evaluar?"
PRINT*,"Para:"
PRINT*,"f(x)=exp(x)-exp(-2x)"
PRINT*,"seleccione 1."
PRINT*,
PRINT*,"Para:"
PRINT*,"f(x)=x-sen(x)"
PRINT*,"seleccione 2."
READ*,p


SELECT CASE(p)

CASE (1)
PRINT*, "Escriba n, el numero de aproximaciones:"
READ*, n
PRINT*, "Escriba x, el numero donde desea evaluar su funcion:"
READ*, x

f_1=0

DO i=1,n

f_1=f_1+(1-(-2)**i)*(x**i)/factorial(i)

END DO

PRINT*, "La funcion evaluada en x, con n aproximaciones es:"
PRINT*, f_1

!-------------------------------------------------------------------

CASE (2)

PRINT*, "Escriba n, el numero de aproximaciones:"
READ*, n
PRINT*, "Escriba x, el numero donde desea evaluar su funcion:"
READ*, x

f_2=0

DO i=1,n

f_2= f_2 + (((-1)**(i+1))*(x**(2*i+1)))/factorial(i)

END DO

PRINT*, "La funcion evaluada en x, con n aproximaciones es:"
PRINT*, f_2

END SELECT


END PROGRAM
