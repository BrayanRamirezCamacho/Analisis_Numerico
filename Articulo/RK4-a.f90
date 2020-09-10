
program RK4

implicit none

real(8),allocatable::x1(:),x2(:),v1(:),v2(:)	!Arreglo que almacena los valores de la función evaluada
real(8),allocatable::t(:)			!Arreglo que almacena los valores del tiempo
real(8)::dx1,dx2,dv1,dv2			!E.D a aproximar
real(8)::ti					!tiempo inicial
real(8)::tf					!tiempo final
real(8)::x1i,x2i,v1i,v2i			!Valor de la función evaluada en el valor inicial
real(8)::h					!Ancho de paso
real(8)::k1,k2,k3,k4,k				!Constructos k que utiliza el método de RK
real(8)::l1,l2,l3,l4,l
real(8)
integer(8)::i,n					!Contadores
integer(8)::error
character(len=80)::m_error	
real(8),parameter::L=14
real(8),parameter::kappa=50
real(8),parameter::m=3
real(8),parameter::omega2=kappa/m

!Le pedimos al usuario el intervalo de aproximación, número de aproximaciones y las C.I.
print*,"Escribe la posición 1 inicial:"
read*, x1i 
print*,"Escribe la posición 2 inicial:"
read*, x2i 
print*,"Escribe la velocidad 1 inicial:"
read*, v1i
print*,"Escribe la velocidad 2 inicial:"
read*, v2i 

ti=0

print*,"Escribe el tiempo final:"
read*,tf
print*,"Escribe el número de puntos, n:"
read*,n

!Alocamos los arreglos según el número de datos:
allocate(t(0:n))
allocate(x1(0:n))
allocate(x2(0:n))
allocate(v1(0:n))
allocate(v2(0:n))

!Establecemos el ancho de paso:
h=(tf-ti)/REAL(n)

!Establecemos los valores inicales:
t(0) = ti
v1(0) = v1i
v2(0) = v2i
x1(0) = x1i
x2(0) = x2i


!Realizamos el algoritmo de RK4 para v2:
do i=1,n
     k1 = dv2(t(i-1) , x1(i-1) , x2(i-1))
     k2 = dv2(t(i-1)+h/2. , x1(i-1)+k1/2. , x2(i-1)+k1/2.)
     k3 = dv2(t(i-1)+h/2. , x1(i-1)+k2/2. , x2(i-1)+k2/2.)
     k4 = dv2(t(i-1)+h , x1(i-1)+k3 , x2(i-1)+k3)
     k = k1 + 2*k2 + 2*k3 + k4
     v2(i) = v2(i-1) + (1./6.)*k*h
     t(i) = t(i-1) + h


end do

!Realizamos el algoritmo de RK4 para x2:
do i=1,n
     k1 = h*dx2(t(i-1),v2(i-1))
     k2 = h*dx2(t(i-1)+h/2.,v2(i-1)+k1/2.)
     k3 = h*dx2(t(i-1)+h/2.,v2(i-1)+k2/2.)
     k4 = h*dx2(t(i-1)+h,v2(i-1)+k3)
     k = k1 + 2*k2 + 2*k3 + k4
     x2(i) = x2(i-1) + (1./6.)*k
     t(i) = t(i-1) + h
end do

!Realizamos el algoritmo de RK4 para v1:
do i=1,n
     k1 = h*dv1(t(i-1),x2(i-1))
     k2 = h*dv1(t(i-1)+h/2.,x2(i-1)+k1/2.)
     k3 = h*dv1(t(i-1)+h/2.,x2(i-1)+k2/2.)
     k4 = h*dv1(t(i-1)+h,x2(i-1)+k3)
     k = k1 + 2*k2 + 2*k3 + k4
     v1(i) = v1(i-1) + (1./6.)*k
     t(i) = t(i-1) + h
end do

!Realizamos el algoritmo de RK4 para x1:
do i=1,n
     k1 = h*dx1(t(i-1),v1(i-1))
     k2 = h*dx1(t(i-1)+h/2.,v1(i-1)+k1/2.)
     k3 = h*dx1(t(i-1)+h/2.,v1(i-1)+k2/2.)
     k4 = h*dx1(t(i-1)+h,v1(i-1)+k3)
     k = k1 + 2*k2 + 2*k3 + k4
     x1(i) = x1(i-1) + (1./6.)*k
     t(i) = t(i-1) + h
end do


!Se generan los archivos de texto con los resultados

print*,"Escriba el nombre del archivo:"
read*,archivo

open( unit=1 , file=archivo , status="new", action="write", iostat=error, iomsg=m_error  )

do i=0,n
	write(1,*)t(i),",",x1(i),",",x2(i),",",v1(i),",",v2(i)
end do

 close(1)

end program RK4

!----------------------------------------------------------------------
!Establecemos las E.D a resolver:
subroutine dv2e(x1,x2)
implicit none
real(8)::L,omega2		!!Argumentos de la ecuación.
real(8),allocatable::x1,x2 
real(8)::dv2			!!E.D a resolver.

  dv2 = omega2 * ( 2*x2 - x1 - L )

end subroutine dv2e

subroutine dx2e(t,v2)
implicit none
real(8),allocatable::t,v2		!!Argumentos de la ecuación.
real(8)::dx2				!!E.D a resolver.

  dx2 = v2 

end subroutine dx2e

subroutine dv1e(t,x2)
implicit none
real(8)::L,omega2		!!Argumentos de la ecuación.
real(8),allocatable::t,x2
real(8)::dv1			!!E.D a resolver.

  dv1 = omega2 * ( L - x2 )

end subroutine dv1e

subroutine dx1e(t,v1)
implicit none
real(8),allocatable::t,v1		!!Argumentos de la ecuación.
real(8)::dx1				!!E.D a resolver.

  dx1 = v1 

end subroutine dx1e











