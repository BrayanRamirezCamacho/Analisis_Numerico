set encoding iso_8859_1
set xlabel "Tiempo (s)"
set ylabel "Posici\363n (m)"
m="./a1"
set title "Posiciones \nrespecto al tiempo"
plot m u 1:2 pt 7 title "x1", m u 1:3 pt 7 title "x2"
