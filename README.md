`vanilla_ipynb`: file that implements the algorithm from "Sampling and homology via bottlenecks" (https://arxiv.org/abs/2010.07976) to generate a point sample and then an updated distance matrix for the projective variety `x[1]*(x[2]^2+x[3]^2+x[4]^2-x[5]^2)-(x[2]^3+x[3]^3+x[4]^3-(1/2)*x[5]^3)` in RP^4 (viewed as a variety in R^5 intersected with the unit sphere).
The distance matrix has roughly size 40000x40000 and about 25GB and is saved as `vanilla_modified_distances.csv`. Store it in the same folder as the file `put_vanilla_modified_distances.csv_here`.

`PH_from_distance_matrix.py` computes the persistent homology of the Vietoris-Rips complex of the distance matrix.
The resulting persistence diagram is saved as `vanilla_plot.png`:

![Persistence diagram for the cubic described above](/PH_from_distance_matrix/vanilla_plot.png)

It has four H1-features, expected is one. 
(There are four dots so close together that it looks like just two dots.
But even two dots would be more than expected.)
It has zero H0-features, expected is one, and it must have at least one.
This suggests that the persistent homology computation is not done correctly.