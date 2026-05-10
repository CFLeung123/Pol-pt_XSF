
---

# Portability Fixes for `parameters.f90`

The original `parameters.f90` relied on several non‑standard extensions of the Intel Fortran compiler that are rejected by gfortran.  
A latent undefined‑behaviour bug was also discovered and corrected.  
This document lists exactly what has been changed in the module.

---

## Changes at a glance

| Issue | Original code (Intel) | Fixed code (gfortran‑compatible) |
|-------|-----------------------|----------------------------------|
| Logical operator `.eqv.` used with integers (`deriv .eqv. .False.` and `deriv .eqv. 0`) | `if(deriv.eqv..False.)then` <br/> `if(deriv.eqv.0)then` | `if(deriv .eq. 0)then` (both branches) |
| Logical value assigned to integer variables (`c = (s0.eq.0)` etc.) | `c = (s0.eq.0)` <br/> `d = (s0.eq.l)` <br/> `e = (s0.eq.(l-1))` | `c = merge(-1, 0, s0 .eq. 0)` <br/> `d = merge(-1, 0, s0 .eq. l)` <br/> `e = merge(-1, 0, s0 .eq. (l-1))` |
| `m_period` function returned an undefined value for `n == 0` | `if (n.ne.0) then` <br/> `m_period = l-n` <br/> `else` <br/> `endif` (no assignment) | `if (n == 0) then` <br/> `m_period = 0` <br/> `else` <br/> `m_period = l - n` <br/> `endif` |

All other subroutines and functions in the module remain untouched.

---

## Detailed description

### 1. Removal of `.eqv.` with integer operands (lines ~140 and ~147)

The variable `deriv` is declared as `integer` (see `input.f90`).  
The Fortran standard requires both operands of `.eqv.` to be `logical`. Intel compilers accept the extension, but gfortran does not.  
The original intent was to test whether `deriv` is zero. Replacing both occurrences with `deriv == 0` preserves the logic and is fully standard compliant.

### 2. Explicit logical‑to‑integer conversion with `merge` (lines ~182‑184)

The original loop‑limits subroutine used

```fortran
c = (s0.eq.0)
d = (s0.eq.l)
e = (s0.eq.(l-1))
```

Intel compilers map `.TRUE.` to `-1` and `.FALSE.` to `0`.  
gfortran rejects this implicit conversion.  
The fix employs the standard intrinsic `merge(-1, 0, logical_condition)`, which guarantees the **same** numerical values (`-1` / `0`) as the Intel behaviour while being accepted by any standard‑compliant compiler.

### 3. Undefined return value in `m_period`

The function `m_period` provides periodic boundary conditions for momentum vectors.  
The original code only assigned a value when `n /= 0`; for `n == 0` the return value was **undefined**.  
Intel compilers often accidentally returned `0`, but gfortran (and many other compilers) may produce arbitrary results.  
This bug caused large discrepancies in several one‑loop diagrams (most notably `_9a` and `_9b`).  

The fix explicitly sets `m_period = 0` when `n == 0`, making the behaviour deterministic and correct across all compilers.

---

## Testing and validation

The correctness of the ported code has been verified by comparing the full output of
gfortran against the original Intel (ifx) results for a small lattice (\(L = 6\),
`mom_deg = 0`, `Gk2_4f` enabled, tree + one‑loop + counterterms).

* **Tree‑level** (`_0`) results are **exactly identical** between the two compilers.
* **Imaginary parts** of every individual one‑loop diagram agree to **at least 14
  significant digits**; most diagrams match to the last printed digit.
* The **total one‑loop** imaginary part (`_1loop`) agrees within \(2\times10^{-15}\).
* **Real parts** of the one‑loop diagrams are theoretically zero; both compilers give
  values of order \(10^{-18}\)–\(10^{-20}\) that differ only because of the order of
  floating‑point additions in OpenMP reductions and minor differences in mathematical
  libraries.  These differences are harmless numerical noise.

The observed deviations are therefore entirely attributable to normal floating‑point
effects.

---

## Compiler flags required for gfortran

These fixes make `parameters.f90` standard Fortran 90/95.  
To compile the whole project, the Makefile should use:

```makefile
f90comp = gfortran
switch  = -O2 -fopenmp -ffree-line-length-0
libs    = -llapack -lblas
```

(`-ffree-line-length-0` is unrelated to this file but is needed for other source files containing long lines.)

---

This document covers only the changes inside `parameters.f90`.  
For the complete build procedure and run instructions, see the main `README.md`.
