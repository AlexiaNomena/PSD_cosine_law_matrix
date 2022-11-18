## PSD cosine law matrix
This repositorie of codes accompanies our paper on: "Qualitative Euclidean Embedding of Disjoint Sets of Points" (link comming soon)

The codes can be used to test the sufficient conditions, provided in our paper above, for obtaining a positive semi-definite matrix **M** with components

**M**<sub>ij</sub> = (*h*(i, 0)<sup>2</sup> + *h*(j, 0)<sup>2</sup>  - *h*(i, j)<sup>2</sup> )/2, where *h* is a specific non-negative symmetric function with defined properties. We term **M** as cosine law matrix because it is obtained from the cosine law between points in Euclidean space.

This is important because **M** is Gram matrix of a set of Euclidean points iff it is positive semi-definite (**M** = **E** **E<sup>T</sup>**). That is any arbitrary set of points generating the cosine law matrix **M** can be embedded into the set of Euclidean points with coordinates obtained on the rows of **E** [(Schoenberg criterion)](https://en.wikipedia.org/wiki/Euclidean_distance_matrix#cite_note-4).

## Description
`working_points.ipynb`: jupyter notebook to demonstrate that the theory works for two sets of points with point coordinated constructed from a contingency dataset and using their representations in Correspondence Analysis.

`rand_points.ipynb`: jupyter notebook to demonstrate that the theory works on randomly generated points (at least it works for the few trials done with our computers)

`rand_points_chgdist.ipynb`: same as `rand_points.ipynb` but the distance within set was multiplied by a positive constant

## Theory
Mathematical statements and proofs are provided in our paper (comming soon)

## License
November 14, 2022: Creative Commons [Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
