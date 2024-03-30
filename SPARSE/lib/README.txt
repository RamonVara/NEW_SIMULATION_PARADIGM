Se ha utilizado la version 5.1 de SUPERLU para la construcción de la 
libreria incluida en EcosimPro, la librería "libsuperlu_5.1.a". 

En el subdirectorio sources_superlu, se proporcionan tres ficheros 
para construir la librería "libsuperlu_5.1.a":

1) d_ecosim_superlu.c: Código fuente de la función de interfaz con EcosimPro
   Este fichero añadir al sub-directorio FORTRAN de superlu. 

2) Makefile para construir la libreria libsuperlu_5.1.a incluyendo el 
   fichero de interfaz con EcosimPro. Este fichero va en el subdirectorio
   FORTRAN de la organización en carpetas de superlu.
   
3) Fichero ..\make.inc que se incluye en el FORTRAN\Makefile. Ha sido preciso 
   cambiar dicho fichero, porque en la compilación estándar de superlu, se obtienen
   2 librerías: superlu y blas. En EcosimPro es preferible juntar estas dos librerías
   en una única librería.
   
El Makefile proporcionado permite construir solo la versión GNU de la librería.





