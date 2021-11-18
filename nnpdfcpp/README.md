# nnpdfcpp
The main programs used in the NNPDF fitting framework.

## Project summary and aim

This project contains the principle fitting code of the NNPDF collaboration,
along with it's associated utility programs.

The main projects are
- filter: Performs fit setup and data kinematic cuts according to an input
  runcard.
- nnfit: Performs a NNPDF replica fit based upon the results directory output
  but filter.

and the extra optional programs:
- evolvefit: evolve PDF neural network using APFEL
- fiatlux: generate a T0 set with the photon PDF via `libfiatlux`
