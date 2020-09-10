!------------------------------------------------------------------------
!Este programa aplica el metodo de Runge-Kutta de orden 4 para
!aproximar numericamente la solucion de un sistema de EDO's acopladas 
!------------------------------------------------------------------------

program RK4

!Programa principal o “manejador”
!Asigna valores para

real::yi !valores iniciales de n variables dependientes
real::ti !valor inicial de la variable independiente
real::tf !valor final de la variable independiente
real::dx !cálculo del tamaño de paso
real::xout !intervalo de salida
real::x
integer::m
real::xpm

read*,y1,y2
xi=0
xf=3600



x = xi
m = 0
xpm = x

do i = 1, 2
	ypim = yi
	yi = yii
end do

do
 xend = x + xout
 
if (xend > xf) then xend = xf
 h = dx
 call Integrator (x, y, 2, h, xend)
 m = m + 1
 xpm = x

 do i = 1, 2
 ypi,m = yi
 end do

 if (x > xf) exit

end do

!OJOJOJOOJJOO
open(unit=1, status="new", action="write", iostat= )

end program

!------------------------------------------------------
!Rutina para tomar un paso de salida
subroutine Integrator (x, y, 2, h, xend)
 do
 if (xend – x < h) then h = xend – x
 call RK4 (x, y, 2, h)
 if (x >= xend) exit

 end do

end subroutine

!------------------------------------------------------
!Método RK de cuarto orden para un sistema de EDO
subroutine RK4(x, y, 4, h)
 call derivadas (x, y, k1)
 
 do i = 1, 4
 ymi = yi + k1i * h/2
 end do
 
 call derivadas (x + h / 2, ym, k2)
 
 do i = 1, 4
 ymi = yi + k2i * h / 2
 end do
 
 call derivadas (x + h / 2, ym, k3)
 
 do i = 1, 4
 yei = yi + k3i * h
 end do

 call derivadas (x + h, ye, k4)
 
 do i = 1, 4
 slopei = (k1i + 2*(k2i+k3i)+k4i)/6
 yi = yi + slopei * h
 end do

 x = x + h

end subroutine

!------------------------------------------------------
!Subrutina para determinar derivadas
subroutine derivadas(x,y,dy)
dy1=
dy2=

end subroutine
