program RK4

implicit none   
integer,parameter::dp = selected_real_kind(32,300)		!Kind
integer,parameter::m = 4 					!número de EDO's acopladas
real(dp)::Y(m)							!Vector para las ecuaciones
real(dp)::ti, tf, t						!Tiempos
real(dp)::h							!Ancho de paso
integer::N, i							!Contadores
real(dp)::w2							!Frecuencia de oscilación
real(dp),parameter::L=14, kappa=50, masa=3			!Constantes
real(dp)::x1i,x2i,v1i,v2i					!Condiciones iniciales
character(len=80)::archivo, m_error,m_error2,m_error3,m_error4	!Nombre del archivo y mensaje de error
integer(dp)::error,error2,error3,error_4			!Errores de lectura
real(dp)::K1,U1,E1,K2,U2,E2					!Energías
real(dp)::P1,P2							!Momentos


!Omega²=k/m
w2 = kappa/masa	 

!Numero de pasos
print*,"Escriba N"
read*,N

!Tiempo inicial
ti = 0

!Tiempo final
print*,"Escriba el tiempo de la simulación:"
read*,tf

!Ancho de paso
h = (tf-ti)/N

!Se pide al usuario las C.I.
print*,"Escriba x1i:"
read*,x1i
print*,"Escriba x2i:"
read*,x2i
print*,"Escriba v1i:"
read*,v1i
print*,"Escriba v2i:"
read*,v2i


!Condiciones iniciales
Y(2) = x1i	! x1(0)
Y(1) = v1i	! v1(0)
Y(4) = x2i	! x2(0)
Y(3) = v2i	! v2(0)


!Se generan los archivos de texto con las posiciones y velocidades

print*,"Escriba el nombre del archivo:"
read*,archivo

open( unit=1 , file=archivo, status="new", action="write", iostat=error, iomsg=m_error  )


!open ( unit=2 , file="energias1.dat" , status="new" , action="write" , iostat=error2 , iomsg=m_error2 )

!open ( unit=3 , file="energías2.dat" , status="new" , action="write" , iostat=error3 , iomsg=m_error3 )

!open ( unit=4 , file="momentos.dat" , status="new" , action="write" , iostat=error4 , iomsg=m_error4 )


t=ti

!Do para calcular los valores
do i = 1,N

!Restricción para las posiciones (específico del problema, resortes de 14m de longitud)
if ( Y(2)>-14.and.Y(2)<28.and.Y(4)<14.and.Y(4)>-28) then

	Y = iterate(t,Y)
	t = t + h

!Energía cinética, potencial y total para 1
!	K1=0.5*masa*Y(1)**2
!	U1=0.5*kappa*Y(2)**2
!	E1 = K1 + U1

!Energía cinética, potencial y total para 2
!	K2=0.5*masa*Y(3)**2
!	U2=0.5*kappa*Y(4)**2
!	E2 = K2 + U2

!Momentos para 1 y 2
!	P1=masa*Y(1)
!	P2=masa*Y(3)


	write(1,*)t,Y(2),Y(4),Y(1),Y(2)

!	write(2,*)t,K1,U1,E1
	
!	write(3,*)t,K2,U2,E2

!	write(4,*)t,P1,P2

else
	exit

end if


end do

!print*, Y

 close(1)

! close(2)

! close(3)

! close(4)


!Graficar en gnuplot, para que se grafique automaticamente, escriba "a1" como nombre del archivo con los datos
call system("gnuplot -p graf.plt")



!-------------------------------------------------
contains

!Función f para calcular el vector f, con la forma de cada ecuacion diferencial

function f(t, Yvec) result (fvec)
real(dp):: t
real(dp):: Yvec(m), fvec(m)

!Yvec(1)=v1
!Yvec(2)=x1
!Yvec(3)=v2
!Yvec(4)=x2

fvec(1) = w2*( L - Yvec(4) ) 		!w2(L-x2) 
fvec(2) = Yvec(1) 			!v1 
fvec(3) = w2*( 2*Yvec(4) - Yvec(2) - L )!w2(2*x2-x1-L)
fvec(4) = Yvec(3) 			!v2

end function

!Función para calcular Y(t_n+1)

function iterate(t, Y_n) result (Y_nplus1)
        real(dp) :: t
        real(dp) :: Y_n(m), Y_nplus1(m)
        real(dp) :: k1(m), k2(m), k3(m), k4(m)

        k1 = h*f(t, Y_n)
        k2 = h*f(t + h/2 , Y_n + k1/2)
        k3 = h*f(t + h/2 , Y_n + k2/2)
        k4 = h*f(t + h , Y_n + k3)

        Y_nplus1 = Y_n + (k1 + 2*k2 + 2*k3 + k4)/6.

end function


end program
