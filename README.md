# rdg-and-fr
A C++ implementation of the robust discontinuous Galerkin (DG) method and the flux reconstruction schemes. This code focuses on the solution of compressible Euler equations and compressible Navier-Stokes equations in 3D space on unstructured hexahedral meshes, although 1D and 2D problems are also included. The parallelization is done by MPI + OpenMP + CUDA.

Code is under construction.

## TODO
- Curvilinear mapping to support curvilinear elements

## Main References

**A.R. Winters, D.A. Kopriva, G.J. Gassner, and F. Hindenlang, "Construction of Modern Robust Nodal Discontinuous Galerkin Spectral Element Methods for the Compressible Navier-Stokes Equations", *Efficient High-Order Discretizations for Computational Fluid Dynamics*, Springer, 2021, pp. 117-196**

**G.J. Gassner, A.R. Winters, and D.A. Kopriva, "Split Form Nodal Discontinuous Galerkin Schemes with Summation-By-Parts Property for the Compressible Euler Equations", *Journal of Computational Physics*, Vol. 327, 2016, pp. 39-66**

**H. Ranocha, M. Schlottke-Lakemper, J. Chan, A.M. Rueda-Ramirez, A.R. Winters, F. Hindenlang, and G.J. Gassner, "Efficient Implementation of Modern Entropy Stable and Kinetic Energy Preserving Discontinuous Galerkin Methods for Conservation Laws", *arXiv:2112.10517[cs.MS]*, https://doi.org/10.48550/arXiv.2112.10517**

**D.A. Kopriva, "Metric Identities and the Discontinuous Spectral Element Method on Curvilinear Meshes", *Journal of Scientific Computing*, Vol. 26, 2006, pp. 301-327**

**H.T. Huynh, "A Flux Reconstruction Approach to High-Order Schemes Including Discontinuous Galerkin Methods", *18th AIAA Computational Fluid Dynamics Conference*, June 2007, Miami, FL, AIAA 2007-4079**

**A. Cicchino, D.C. Del Rey Fernandez, S. Nadarajah, J. Chan, and M.H. Carpenter, "Provably Stable Flux Reconstruction High-Order Methods on Curvilinear Elements", *Journal of Computational Physics*, Vol. 463, 2022, 111259**
