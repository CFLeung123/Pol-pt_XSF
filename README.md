Below is a README.md file that explains how to build and run the code with gfortran on a Linux-like system (e.g., WSL). You can place it in the root directory of the project.

```markdown
# pt_XSF

Fortran 90/95 code for the perturbative computation of correlation functions
in the chirally rotated Schrödinger functional (ΧSF).  
Originally written for the Intel compiler, it can be easily adapted to
GNU Fortran (gfortran).

## Overview

The program computes tree‑level, one‑loop, and counterterm contributions to
several two‑ and four‑fermion correlation functions on a finite lattice.
It uses:

* Wilson/clover fermions with ΧSF boundary conditions,
* plaquette or Lüscher–Weisz gauge action,
* analytic Feynman rules (vertex functions) and numerical inversion of
  the Dirac operator (LAPACK banded matrix routines).

Parallelisation is achieved with OpenMP (`!$OMP PARALLEL DO`) for the
most expensive momentum loops.

## Requirements

* **gfortran** (or any Fortran 90/95 compiler; Intel ifort works as well)
* **LAPACK** and **BLAS** libraries (e.g., `liblapack-dev libblas-dev`)
* **OpenMP** support (usually included in gfortran)
* **Bash** shell for the helper scripts

## Directory structure

```
pt_XSF/
├── source/             # Fortran sources + Makefile
├── in.dat.src          # Input file template
├── sizelist            # List of lattice sizes for batch runs (optional)
├── compile.sh          # Shell script to compile the program
├── pt_runner.sh        # Shell script for batch production runs
├── bin/                # Created by compile.sh – executable + mod files
├── output/             # Created by pt_runner.sh – numerical results
└── README.md           # This file
```

## Building

### 1. Install prerequisites (Ubuntu / WSL)

```bash
sudo apt update
sudo apt install gfortran liblapack-dev libblas-dev
```

### 2. Adapt the source for gfortran

The code uses two non‑standard Fortran idioms that are accepted by Intel
compilers but rejected by gfortran. Before compiling, apply the following
changes to `source/parameters.f90`.

**(a)** Replace logical‑integer equality (around line 140 and 147):

```fortran
! OLD
if (deriv .eqv. .False.) then   ! line ~140
if (deriv .eqv. 0) then         ! line ~147

! NEW
if (deriv == 0) then
if (deriv == 0) then
```

**(b)** Replace logical‑to‑integer conversions (lines ~182‑184) to avoid
warnings (optional but clean). For example:

```fortran
! OLD
c = (s0.eq.0)
d = (s0.eq.l)
e = (s0.eq.(l-1))

! NEW
c = 0; d = 0; e = 0
if (s0 == 0)   c = 1
if (s0 == l)   d = 1
if (s0 == l-1) e = 1
```

### 3. Modify the Makefile

Open `source/Makefile` and set:

```makefile
f90comp = gfortran
switch  = -O2 -fopenmp -ffree-line-length-0 -Wno-conversion
libs    = -llapack -lblas
```

The flag `-ffree-line-length-0` is required because some statements are
longer than 132 characters.

### 4. Update the helper scripts

Both `compile.sh` and `pt_runner.sh` contain a line that tries to load
Intel compiler variables. **Remove or comment out** that line:

```bash
# source /opt/intel/bin/compilervars.sh intel64
```

(If you are using ifort, keep the line – but the rest of this README
assumes gfortran.)

If you intend to use more than one CPU core, also adjust the export in
`pt_runner.sh`:

```bash
export OMP_NUM_THREADS=12      # example for a 14‑core machine
```

### 5. Compile

From the project root directory:

```bash
cd source && make clean && cd ..   # start clean
bash compile.sh                    # build the executable
```

After success, the executable `pt_simple.out` and the Fortran module
files (`*.mod`) will be inside `bin/`.

## Preparing the input file

The file `in.dat.src` is a template. The keyword `=CFG=` will be replaced
by the actual lattice size.

* **For a quick test:** manually create `bin/in.dat`
  ```bash
  sed 's/=CFG=/6/' in.dat.src > bin/in.dat
  ```

* A minimal `in.dat` that calculates only the tree‑level `g1` boundary
  correlator (very fast) looks like:
  ```
  Setup_parameters:
  ****************
  l_size  6
  m0      0.0
  …       (other parameters kept default)
  mom_deg 1
  …
  What to compute:
  ****************
  g1      1
  …
  tree    1
  1loop   0
  counter 0
  ```

  Adjust the switches for your physics case.

## Running

### Single run (foreground, output to screen)

```bash
cd bin
./pt_simple.out
```

### Batch run with multiple sizes

1. Create a file `sizelist` containing one lattice size per line, e.g.:
   ```
   6
   8
   12
   ```
2. Edit `pt_runner.sh` – make sure the Intel line is removed and
   `OMP_NUM_THREADS` is set appropriately.
3. Run:
   ```bash
   bash pt_runner.sh
   ```
   For each size, the script replaces `=CFG=` in `in.dat.src`,
   executes `pt_simple.out`, and saves the output as
   `output/<size>.dat`.

## Output

The program prints the correlation functions in plain text. Each line
contains:

* observable label (e.g., `gV0`, `lA_`, `Gk2_`)
* flavour combination (`uu`, `ud`, …)
* diagram label (`_0` tree, `_1a` … `_13b` one‑loop diagrams,
  `_1loop` total one‑loop, `_dm`, `_zf`, `_ds` counterterms)
* time‑slice index `i` (only for time‑dependent observables)
* real and imaginary parts of the complex result

Example:

```
gV0 ud _0 : 5  1.2345678E-02  -2.3456789E-03
Gk2_ II____ g5____ uu ud _1a= 3  -1.234E-03   2.345E-04
```

## Notes on performance

* The most expensive diagrams are the four‑fermion correlation functions
  (`G1_4f`, `G2_4f`, `Gk1_4f`, `Gk2_4f`). Their computation scales
  roughly as `O(L^9)` without momentum degeneracy.
* To obtain the fastest output, enable `mom_deg = 1` in the input file.
  It exploits momentum degeneracy and reduces the number of loop
  configurations by about a factor of 6.
* The OpenMP parallelisation benefits multi‑core systems. Set
  `OMP_NUM_THREADS` to a value slightly less than the number of physical
  cores for best performance.

## Original environment

The code was originally developed for the Intel Fortran compiler (`ifort`)
with `source /opt/intel/bin/compilervars.sh intel64`. All non‑standard
extensions have been removed (or are trivial to adapt) so that the program
compiles cleanly with gfortran.
```

This README summarises everything needed to compile and run the code using gfortran, based on our step‑by‑step debugging.
