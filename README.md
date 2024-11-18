# Code accompanying the manuscript "Global stability of ecological and evolutionary dynamics via equivalence", by Stefano Allesina

The repository contains code that can find weights for the Lyapunov function analyzed in the manuscript. In particular:

- `RE_n_3.nb` is a Mathematica Notebook containing code for finding weights for an arbitrary 3-dimensional Replicator Equation with a globally stable equilibrium in the interior of the simplex. A `pdf` version is also provided.
- The file `find_weights.R` contains code that searches numerically for appropriate weights of a matrix of arbitrary size. It uses the file `check_copositivity.R` that contains code for testing whether a matrix is copositive.

All the examples presented in the manuscript are also analyzed in the code presented here. 

For any question/comment/clarification, contact Stefano Allesina `sallesina@uchicago.edu`
