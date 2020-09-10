  !!El objetivo de este programa es el de aproximar la solución de una
  !!ecuación diferencial (E.D) utilizado el método de Runge-Kutta de
  !!orden 4 (RK4).
  !!-------------------------------------------------------------------

PROGRAM rungekuta4

  IMPLICIT NONE
  !Diccionario de variables
  REAL(8),ALLOCATABLE::y(:)     !!Arreglo que almacena los valores de la función evaluada.
  REAL(8),ALLOCATABLE::t(:)     !!Arreglo que almacena los valores a evaluar.
  REAL(8)::ed                   !!E.D a la cual aproximaremos su resultado.
  REAL(8)::ti                   !!Valor inicial.
  REAL(8)::yi                   !!Valor de la función evaluada en el valor inicial.
  REAL(8)::h                    !!Ancho de paso.
  REAL(8)::k1,k2,k3,k4          !!Constructos k que utiliza el método de RK.
  REAL(8)::k                    !!Sumatoria de las k.
  REAL(8)::a,b                  !!Intervalo de la E.D.
  INTEGER::i,n                  !!Contadores.

  !Explicamos al usuario las funciones del programa:
  PRINT*,"Bienvenido al programa de Runge-Kutta, donde utilizaremos el método"
  PRINT*,"de Runge-Kutta de grado 4 (RK4) para resolver una ecuación diferencial."
  PRINT*,"En este caso, la E.D a resolver es:"
  PRINT*,"y' = 1 + y**2 + t**3"

    !Le pedimos al usuario el intervalo de aproximación, número de aproximaciones y las C.I.
  PRINT*,"Introduce el número de datos que deseas aproximar:"
  READ*,n
  
     PRINT*,"Introduce el inicio del intervalo de aproximación:"
     READ*,a
     PRINT*,"Introduce el final del intervalo de aproximación:"
     READ*,b
     IF(a>b)THEN
        PRINT*,"Lo sentimos, pero el valor inicial no puede ser mayor que el valor final."
        PRINT*,"Por favor, intentalo de nuevo."
        CONTINUE
     ELSE
        EXIT
     END IF
  
  PRINT*,"Introduce el valor inicial a evaluar:"
  READ*,ti
  PRINT*,"Introduce el valor de la función evaluada en el valor inicial:"
  READ*,yi

  !Alocamos los arreglos según el número de datos:
  ALLOCATE(t(0:n))
  ALLOCATE(y(0:n))

  !Establecemos el ancho de paso:
  h=(b-a)/REAL(n)
  
  !Establecemos los valores inicales:
  y(0) = yi
  t(0) = ti

    
  !Realizamos el algoritmo de RK4:
  DO i=1,n
     k1 = h*ed(t(i-1),y(i-1))
     k2 = h*ed(t(i-1)+h/2.,y(i-1)+k1/2.)
     k3 = h*ed(t(i-1)+h/2.,y(i-1)+k2/2.)
     k4 = h*ed(t(i-1)+h,y(i-1)+k3)
     k = k1 + 2*k2 + 2*k3 + k4
     y(i) = y(i-1) + (1./6.)*k
     t(i) = t(i-1) + h
  END DO
  !En el algoritmo anterior, los constructos K son, geométricamente,
  !rectas que utilizamos para aproximar la siguiente iteración.

    !Mostramos al usuario los resultados:"
  PRINT*,"Tu tabla de datos resultante es:"
  DO i=0,n
     PRINT*,t(i),y(i)
  END DO
  
  
END PROGRAM rungekuta4

!Establecemos la E.D a resolver:
FUNCTION ed(t,y)
  REAL(8)::t,y                  !!Argumentos de la ecuación.
  REAL(8)::ed                   !!E.D a resolver.

  ed = 1 + y**2 + t**3

END FUNCTION ed
