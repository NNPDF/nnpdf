# Database of theories implemented in the NNPDF framework

## General guidelines

Each unique theory is assigned an integer as an identifier of the form `XX_YYY_ZZZ` where

- `XX` refers to the timeframe of the project under which the theory is created/used, for this we utilize the releases of unpolarized global PDF fits.
- `YYY` refers to the project or family of theories. For instance, `000` is intrinsic charm, `001` perturbative charm, etc.
- `ZZZ` is an index with no specific meaning for the theories within a project. In the unlikely case that a project needs more than 1000 theories, two different `YYY` can refer to the same project.

It is common for different projects to share a subset of theories. In those cases, links may be used to save space (specially in the server), although they would be distinct theories when downloading them to a local computer or a cluster.

## Theories in use

This is a (non-exhaustive) index of theories currently implemented.

### 40

- `40_000_000`: NNPDF4.0 baseline
- `_000_`: default 4.0 theories, nFONLL with intrinsic charm.
- `_001_`: 4.0 theories, nFONLL with perturbative charm.
- `_002_`: 4.0 theories, nFONLL, N3LO

### 41
- `41_000_000`: NNPDF4.1 baseline
- `_100_`: Theories for NNPDFpol2.0

## Old theories
Theories below <= 1026 correspond to fits for the NNPDF4.0 global analysis (and before)
