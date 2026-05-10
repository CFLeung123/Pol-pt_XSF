Below is a **README** focused on the changes made to `parameters.f90` to achieve full portability between Intel (`ifort`/`ifx`) and GNU (`gfortran`) compilers, along with a fix for an undefined‑behaviour bug.

---

# Portability Fixes for `parameters.f90`

The original `parameters.f90` relied on several non‑standard extensions of the Intel Fortran compiler that are not accepted by gfortran.  
Additionally, a subtle uninitialised variable bug was discovered during the cross‑compiler validation.  
All modifications are collected in this document.

## Changes at a glance

| Issue | Original code | Fixed code |
|-------|---------------|------------|
| Logical operator `.eqv.` used with integers (`deriv.eqv..False.`) | `if (deriv.eqv..False.) then` <br/> `if (deriv.eqv.0) then` | `if (deriv .eq. 0) then` (both branches) |
| Logical value assigned to integer (`c = (s0.eq.0)`) | `c = (s0.eq.0)` <br/> `d = (s0.eq.l)` <br/> `e = (s0.eq.(l-1))` | `c = merge(-1, 0, s0 .eq. 0)` <br/> `d = merge(-1, 0, s0 .eq. l)` <br/> `e = merge(-1, 0, s0 .eq. (l-1))` |
| `even_or_odd` function used floating point instead of modulo arithmetic | `a = (-1.0)**l` <br/> `if(a.gt.0.0) then` | `if (mod(l, 2) == 0) then` |
| `m_period` function returned an undefined value for `n == 0` | `if(n.ne.0) then` <br/> `m_period = l-n` <br/> `else` <br/> `endif` (no assignment in else branch) | `if (n == 0) then` <br/> `m_period = 0` <br/> `else` <br/> `m_period = l - n` <br/> `endif` |

All other subroutines and functions remain unchanged.

## Detailed explanation

### 1. Removal of `.eqv.` with integer operands (line ~140, ~147)

The variable `deriv` is declared as `integer` in `input.f90`.  
The Fortran standard requires both operands of `.eqv.` to be of type `logical`.  
Intel compilers accept the extension, gfortran (by default) does not.  
The original logic was a simple zero/non‑zero test; replacing it with `deriv == 0` preserves the same behaviour and is standard conformant.

### 2. Logical‑to‑integer conversion in `loop_limits` (lines ~182‑184)

The original code used:

```fortran
c = (s0.eq.0)
d = (s0.eq.l)
e = (s0.eq.(l-1))
```

Intel compilers map `.TRUE.` to `-1` and `.FALSE.` to `0`.  
gfortran rejects this implicit conversion.  
The fix uses the standard intrinsic `merge(-1, 0, logical_condition)`, which explicitly returns `-1` for `.TRUE.` and `0` for `.FALSE.`.  
This guarantees identical arithmetic in the subsequent index‑boundary expressions.

### 3. Safer implementation of `even_or_odd`

The old implementation used `a = (-1.0)**l` and checked `a > 0.0`.  
Although it worked, floating‑point exponentiation is unnecessarily fragile and can lead to warnings with some compiler flags.  
The replacement `if (mod(l, 2) == 0)` is a pure integer operation and is standard conformant across all Fortran compilers.

### 4. Critical fix in `m_period` – undefined return value

The `m_period` function implements periodic boundary conditions for momenta.  
The original code only assigned a value when `n /= 0`; for `n == 0` the function returned whatever happened to be in the memory location (undefined behaviour).  
Intel compilers often returned `0` by accident, but gfortran produced a different (sometimes very large) value.  
This caused large discrepancies in certain diagrams (e.g. `_9a`, `_9b`) that involve `mqloop` built from `m_period`.  
The fix explicitly sets `m_period = 0` when `n == 0`.

## How to verify the correctness

After applying these changes and recompiling with gfortran:

* Tree‑level results should be **exactly identical** to those obtained with the Intel compiler.
* One‑loop diagrams should agree to **machine precision** in the imaginary parts.  
  (Small differences of order `1e-15` may remain in the real parts for some diagrams due to different floating‑point reduction order in OpenMP and different math libraries; these are harmless.)
* The total one‑loop result `_1loop` must match within similar tolerance.

## Compilation with gfortran

The modified `parameters.f90` is now fully standard Fortran 90/95.  
Use the following flags in the Makefile to compile without errors:

```makefile
f90comp = gfortran
switch  = -O2 -fopenmp -ffree-line-length-0
libs    = -llapack -lblas
```

`-ffree-line-length-0` is not related to this file, but is necessary for other source files that contain long lines.

---

This file documents only the content of the `parameters.f90` module.  
For a complete build and run guide, refer to the main `README.md`.
