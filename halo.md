J'ai choisi de ne séparer la grille qu'en hauteur:

ny
^
|
|
|
0- - - - > nx


Pour N (>= 1) le nombre de rang mpi:

Si N = 1 :
0 ou Deux halos, 0, nx

Si N = 2:
Quatre halos, 0, nx//2-1, nx//2, nx


Si N = 3 : 
Six halos, 0, nx//3-1, nx//3, 2*(nx//3) - 1, 2*(nx//3)-1 
