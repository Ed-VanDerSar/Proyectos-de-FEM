# Proyectos-de-FEM
Proyectos realizado para el [curso de implementación y teoría del Método del Elemento](https://danielcq-math.github.io/cursos/mfem_2023_II/index.html) Finito del Posgrado en Ciencias Matemáticas de la UNAM.

Se usa el método del elemento finito de Lagarange y Crouzeix–Raviart (de primer orden) aproximar la solución débil al problema de Poisson $-\Nabla\cdot(\kappa \nabla u)=f$ en $\Omega=(0,1)$ y $\Omega=(0,1)\times(0,1)$ con condiciones de forntera de Dirichlet. El algoritmo está documentado. El imput es la función dada $f$ y el output son las trablas de convergencia en norma de Legesgue $L^2$ y de Sobolev $H^1$ con las respectivas tasas de convergencia.
