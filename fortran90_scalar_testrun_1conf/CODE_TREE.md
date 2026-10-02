# Code Tree

Generated from the Fortran sources in this directory. Dependencies below mean direct `CALL` relationships and shared `COMMON` block usage detected from the source.

## Repository Tree

```text
axion/
├── torelons_4D_adjoint.f        # main torelon-correlation analysis code
├── luesher.f                    # Luescher/Wolff random number generator
└── torelons_4D_adjoint          # compiled Mach-O arm64 executable
```

## Build-Level Dependencies

- `torelons_4D_adjoint.f` defines the implicit main program and most numerical routines.
- `luesher.f` provides `RINI`, `RNDNUM`, and `MULTI_RNDNUM`; the main program currently calls `RINI`.
- No `include` files or Fortran `module` dependencies were detected.
- Runtime/library usage includes `cpu_time`, `execute_command_line`, file I/O, and Fortran numeric intrinsics such as `sqrt`, `log`, `exp`, `dreal`, `dimag`, `dcmplx`, `sngl`, and `mod`.

## Program Call Tree

- `MAIN` - implicit program, torelons_4D_adjoint.f:1
  - `SETUP` - subroutine, torelons_4D_adjoint.f:18168
  - `RINI` - subroutine, luesher.f:175
  - `READ_GF` - subroutine, torelons_4D_adjoint.f:251
  - `ACTION` - subroutine, torelons_4D_adjoint.f:17755
    - `VMX` - subroutine, torelons_4D_adjoint.f:18431
    - `HERM` - subroutine, torelons_4D_adjoint.f:18474
    - `TRVMX` - subroutine, torelons_4D_adjoint.f:18454
    - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
  - `POLYLOOP` - subroutine, torelons_4D_adjoint.f:391
    - `VMX` - subroutine, torelons_4D_adjoint.f:18431
    - `TRVMX` - subroutine, torelons_4D_adjoint.f:18454
  - `MEASURE` - subroutine, torelons_4D_adjoint.f:888
    - `SETUPB` - subroutine, torelons_4D_adjoint.f:17678
    - `BLOCK` - subroutine, torelons_4D_adjoint.f:2361
      - `SMEAR1` - subroutine, torelons_4D_adjoint.f:2440
        - `DIAG` - subroutine, torelons_4D_adjoint.f:2619
          - `VMX` - subroutine, torelons_4D_adjoint.f:18431
          - `HERM` - subroutine, torelons_4D_adjoint.f:18474
        - `VMX` - subroutine, torelons_4D_adjoint.f:18431
        - `HERM` - subroutine, torelons_4D_adjoint.f:18474
        - `RENORMBS` - subroutine, torelons_4D_adjoint.f:2708
          - `VMX` - subroutine, torelons_4D_adjoint.f:18431
          - `RUNGRB` - subroutine, torelons_4D_adjoint.f:2779
            - `SUBGRB` - subroutine, torelons_4D_adjoint.f:2805
        - `DODET` - subroutine, torelons_4D_adjoint.f:2570
          - `DETNANT` - subroutine, torelons_4D_adjoint.f:18063
            - `DETMAT` - subroutine, torelons_4D_adjoint.f:18080
              - `LUDCMP` - subroutine, torelons_4D_adjoint.f:18098
      - `VMX` - subroutine, torelons_4D_adjoint.f:18431
      - `THERML1` - subroutine, torelons_4D_adjoint.f:2849
        - `VMX` - subroutine, torelons_4D_adjoint.f:18431
        - `HERM` - subroutine, torelons_4D_adjoint.f:18474
        - `RENORMBS` - subroutine, torelons_4D_adjoint.f:2708
          - `VMX` - subroutine, torelons_4D_adjoint.f:18431
          - `RUNGRB` - subroutine, torelons_4D_adjoint.f:2779
            - `SUBGRB` - subroutine, torelons_4D_adjoint.f:2805
        - `DETNANT` - subroutine, torelons_4D_adjoint.f:18063
          - `DETMAT` - subroutine, torelons_4D_adjoint.f:18080
            - `LUDCMP` - subroutine, torelons_4D_adjoint.f:18098
      - `POT` - subroutine, torelons_4D_adjoint.f:12345
    - `MLINE1` - subroutine, torelons_4D_adjoint.f:922
      - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
      - `JACK` - subroutine, torelons_4D_adjoint.f:18296
        - `FITT` - subroutine, torelons_4D_adjoint.f:18385
    - `MLINE2` - subroutine, torelons_4D_adjoint.f:1043
      - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
      - `JACK` - subroutine, torelons_4D_adjoint.f:18296
        - `FITT` - subroutine, torelons_4D_adjoint.f:18385
    - `MLINE3` - subroutine, torelons_4D_adjoint.f:1169
      - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
      - `JACK` - subroutine, torelons_4D_adjoint.f:18296
        - `FITT` - subroutine, torelons_4D_adjoint.f:18385
    - `MLINE4` - subroutine, torelons_4D_adjoint.f:1675
      - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
      - `JACK` - subroutine, torelons_4D_adjoint.f:18296
        - `FITT` - subroutine, torelons_4D_adjoint.f:18385
    - `MLINE5` - subroutine, torelons_4D_adjoint.f:2018
      - `JACKM` - subroutine, torelons_4D_adjoint.f:18224
      - `JACK` - subroutine, torelons_4D_adjoint.f:18296
        - `FITT` - subroutine, torelons_4D_adjoint.f:18385
  - `TODISK1` - subroutine, torelons_4D_adjoint.f:449
  - `TODISK2` - subroutine, torelons_4D_adjoint.f:697

## Direct Routine Dependencies

### `MAIN`

- Location: `torelons_4D_adjoint.f:1`
- Kind/signature: `implicit main program`
- Description: Drives the full analysis: initializes geometry/random state, reads each gauge configuration, measures observables, and writes accumulated results.
- Calls: `SETUP`, `RINI`, `READ_GF`, `ACTION`, `POLYLOOP`, `MEASURE`, `TODISK1`, `TODISK2`
- Called by: none
- Shared state: `/ARRAYS/`

### `QUATERNION_TO_MATRIX`

- Location: `torelons_4D_adjoint.f:235`
- Kind/signature: `function quaternion_to_matrix(quaternion)`
- Description: Converts a four-component quaternion into a `2 x 2` complex SU(2)-style matrix.
- Calls: none
- Called by: none
- Shared state: none

### `READ_GF`

- Location: `torelons_4D_adjoint.f:251`
- Kind/signature: `SUBROUTINE READ_GF(filename)`
- Description: Reads a stream-format gauge field file, converts quaternion links to complex matrices, and stores them in `/ARRAYS/`.
- Calls: none
- Called by: `MAIN`
- Shared state: `/ARRAYS/`

### `READ_GAUGE_ETMC`

- Location: `torelons_4D_adjoint.f:316`
- Kind/signature: `SUBROUTINE READ_GAUGE_ETMC(filename) ! to be fixed !`
- Description: Alternative ETMC/direct-access gauge reader that loads complex links into `/ARRAYS/`; marked in-source as unfinished.
- Calls: none
- Called by: none
- Shared state: `/ARRAYS/`

### `POLYLOOP`

- Location: `torelons_4D_adjoint.f:391`
- Kind/signature: `SUBROUTINE POLYLOOP()`
- Description: Computes Polyakov-loop-like traces in all four directions from the gauge links and nearest-neighbor tables.
- Calls: `VMX`, `TRVMX`
- Called by: `MAIN`
- Shared state: `/ARRAYS/`, `/NEXT/`

### `TODISK1`

- Location: `torelons_4D_adjoint.f:449`
- Kind/signature: `SUBROUTINE TODISK1`
- Description: Writes zero-momentum/binned correlator, Polyakov-line, plaquette, and action data from shared measurement blocks.
- Calls: none
- Called by: `MAIN`
- Shared state: `/PBLOCKP/`, `/PBLOCKJ0PP/`, `/PBLOCKJ0PM/`, `/PBLOCKJ0MP/`, `/PBLOCKJ0MM/`, `/PBLOCKJ1P/`, `/PBLOCKJ1M/`, `/PBLOCKJ2PP/`, `/PBLOCKJ2PM/`, `/PBLOCKJ2MP/`, `/PBLOCKJ2MM/`, `/PBLOCKJ0/`, `/PBLOCKJ1/`, `/PBLOCKJ2/`, `/PBLOCKJALL/`, `/POLYT/`, `/ASTORE/`

### `TODISK2`

- Location: `torelons_4D_adjoint.f:697`
- Kind/signature: `SUBROUTINE TODISK2`
- Description: Writes nonzero-momentum correlator blocks and derived channel data from the momentum measurement arrays.
- Calls: none
- Called by: `MAIN`
- Shared state: `/PBLOCKMOMP/`, `/PBLOCKMOMJ0P/`, `/PBLOCKMOMJ0M/`, `/PBLOCKMOMJ1/`, `/PBLOCKMOMJ2P/`, `/PBLOCKMOMJ2M/`, `/PBLOCKMOMJ0/`, `/PBLOCKMOMJ2/`, `/PBLOCKMOMJALL/`

### `MEASURE`

- Location: `torelons_4D_adjoint.f:888`
- Kind/signature: `SUBROUTINE MEASURE(ITER,NTOT)`
- Description: Coordinates one measurement pass, managing bin counters and calling blocking plus the five line-analysis routines.
- Calls: `SETUPB`, `BLOCK`, `MLINE1`, `MLINE2`, `MLINE3`, `MLINE4`, `MLINE5`
- Called by: `MAIN`
- Shared state: `/ITEM/`

### `MLINE1`

- Location: `torelons_4D_adjoint.f:922`
- Kind/signature: `SUBROUTINE MLINE1`
- Description: Analyzes the primary zero-momentum Polyakov-line correlator block and estimates averages/errors with jackknife routines.
- Calls: `JACKM`, `JACK`
- Called by: `MEASURE`
- Shared state: `/PBLOCKP/`, `/ITEM/`, `/FIT/`

### `MLINE2`

- Location: `torelons_4D_adjoint.f:1043`
- Kind/signature: `SUBROUTINE MLINE2`
- Description: Analyzes momentum-resolved Polyakov-line correlators and jackknife-derived fit quantities.
- Calls: `JACKM`, `JACK`
- Called by: `MEASURE`
- Shared state: `/PBLOCKMOMP/`, `/ITEM/`, `/FIT/`

### `MLINE3`

- Location: `torelons_4D_adjoint.f:1169`
- Kind/signature: `SUBROUTINE MLINE3`
- Description: Processes zero-momentum correlator channels split by angular-momentum/parity/reflection operator classes.
- Calls: `JACKM`, `JACK`
- Called by: `MEASURE`
- Shared state: `/PBLOCKP/`, `/PBLOCKJ0PP/`, `/PBLOCKJ0PM/`, `/PBLOCKJ0MP/`, `/PBLOCKJ0MM/`, `/PBLOCKJ1P/`, `/PBLOCKJ1M/`, `/PBLOCKJ2PP/`, `/PBLOCKJ2PM/`, `/PBLOCKJ2MP/`, `/PBLOCKJ2MM/`, `/ITEM/`, `/FIT/`

### `MLINE4`

- Location: `torelons_4D_adjoint.f:1675`
- Kind/signature: `SUBROUTINE MLINE4`
- Description: Processes mixed zero- and finite-momentum channel correlators for selected operator classes.
- Calls: `JACKM`, `JACK`
- Called by: `MEASURE`
- Shared state: `/PBLOCKP/`, `/PBLOCKJ0PP/`, `/PBLOCKJ0PM/`, `/PBLOCKJ0MP/`, `/PBLOCKJ0MM/`, `/PBLOCKJ1P/`, `/PBLOCKJ1M/`, `/PBLOCKJ2PP/`, `/PBLOCKJ2PM/`, `/PBLOCKJ2MP/`, `/PBLOCKJ2MM/`, `/PBLOCKMOMJ0P/`, `/PBLOCKMOMJ0M/`, `/PBLOCKMOMJ1/`, `/PBLOCKMOMJ2P/`, `/PBLOCKMOMJ2M/`, `/PBLOCKMOMJ0/`, `/PBLOCKMOMJ2/`, `/ITEM/`, `/FIT/`

### `MLINE5`

- Location: `torelons_4D_adjoint.f:2018`
- Kind/signature: `SUBROUTINE MLINE5`
- Description: Processes the remaining mixed channel correlators, mirroring `MLINE4` for additional momentum/operator sectors.
- Calls: `JACKM`, `JACK`
- Called by: `MEASURE`
- Shared state: `/PBLOCKP/`, `/PBLOCKJ0PP/`, `/PBLOCKJ0PM/`, `/PBLOCKJ0MP/`, `/PBLOCKJ0MM/`, `/PBLOCKJ1P/`, `/PBLOCKJ1M/`, `/PBLOCKJ2PP/`, `/PBLOCKJ2PM/`, `/PBLOCKJ2MP/`, `/PBLOCKJ2MM/`, `/PBLOCKMOMJ0P/`, `/PBLOCKMOMJ0M/`, `/PBLOCKMOMJ1/`, `/PBLOCKMOMJ2P/`, `/PBLOCKMOMJ2M/`, `/PBLOCKMOMJ0/`, `/PBLOCKMOMJ2/`, `/ITEM/`, `/FIT/`

### `BLOCK`

- Location: `torelons_4D_adjoint.f:2361`
- Kind/signature: `SUBROUTINE BLOCK`
- Description: Builds blocked/smeared spatial links, thermal line operators, and potential/correlation accumulators for each time slice.
- Calls: `SMEAR1`, `VMX`, `THERML1`, `POT`
- Called by: `MEASURE`
- Shared state: `/ARRAYS/`, `/ARRAYB/`, `/ASMEAR1/`, `/ASMEAR2/`, `/NEXTB/`

### `SMEAR1`

- Location: `torelons_4D_adjoint.f:2440`
- Kind/signature: `SUBROUTINE SMEAR1`
- Description: Performs one smearing step on spatial links using staples, matrix products, renormalization, and determinant normalization.
- Calls: `DIAG`, `VMX`, `HERM`, `RENORMBS`, `DODET`
- Called by: `BLOCK`
- Shared state: `/ASMEAR1/`, `/ASMEAR2/`, `/DIAGIN/`, `/DIAGOUT/`

### `DODET`

- Location: `torelons_4D_adjoint.f:2570`
- Kind/signature: `SUBROUTINE DODET`
- Description: Normalizes smeared link matrices by their determinant to project them back toward the target matrix group.
- Calls: `DETNANT`
- Called by: `SMEAR1`
- Shared state: `/ASMEAR1/`

### `DIAG`

- Location: `torelons_4D_adjoint.f:2619`
- Kind/signature: `SUBROUTINE DIAG(MU)`
- Description: Builds diagonal/staple-like transport products for direction `MU`, using neighbor shifts and Hermitian conjugates.
- Calls: `VMX`, `HERM`
- Called by: `SMEAR1`
- Shared state: `/DIAGOUT/`, `/DIAGIN/`

### `RENORMBS`

- Location: `torelons_4D_adjoint.f:2708`
- Kind/signature: `SUBROUTINE RENORMBS(UU1,UREN11)`
- Description: Renormalizes a blocked/smeared matrix through Gram-Schmidt-like orthogonalization and matrix reconstruction.
- Calls: `VMX`, `RUNGRB`
- Called by: `SMEAR1`, `THERML1`
- Shared state: none

### `RUNGRB`

- Location: `torelons_4D_adjoint.f:2779`
- Kind/signature: `SUBROUTINE RUNGRB(B11,C11)`
- Description: Runs the row/column sub-block updates used by the blocked-link renormalization step.
- Calls: `SUBGRB`
- Called by: `RENORMBS`
- Shared state: `/SUBB/`

### `SUBGRB`

- Location: `torelons_4D_adjoint.f:2805`
- Kind/signature: `SUBROUTINE SUBGRB(LDU,LDL,B11,C11)`
- Description: Applies one elementary sub-block rotation/update within the Gram-Schmidt-like renormalization routine.
- Calls: none
- Called by: `RUNGRB`
- Shared state: `/SUBB/`

### `THERML1`

- Location: `torelons_4D_adjoint.f:2849`
- Kind/signature: `SUBROUTINE THERML1(N4,IBLL)`
- Description: Constructs thermal/torelon line operators for a time slice and blocking level, including momentum-resolved variants.
- Calls: `VMX`, `HERM`, `RENORMBS`, `DETNANT`
- Called by: `BLOCK`
- Shared state: `/LINES/`, `/LINESMOM/`, `/ARRAYB/`, `/NEXTB/`, `/DUMMY/`

### `POT`

- Location: `torelons_4D_adjoint.f:12345`
- Kind/signature: `SUBROUTINE POT`
- Description: Accumulates correlation matrices from line operators across operator channels, momenta, bins, and time separations.
- Calls: none
- Called by: `BLOCK`
- Shared state: `/PBLOCKP/`, `/PBLOCKMOMP/`, `/PBLOCKJ0PP/`, `/PBLOCKJ0PM/`, `/PBLOCKJ0MP/`, `/PBLOCKJ0MM/`, `/PBLOCKJ1P/`, `/PBLOCKJ1M/`, `/PBLOCKJ2PP/`, `/PBLOCKJ2PM/`, `/PBLOCKJ2MP/`, `/PBLOCKJ2MM/`, `/PBLOCKJ0/`, `/PBLOCKJ1/`, `/PBLOCKJ2/`, `/PBLOCKJALL/`, `/PBLOCKMOMJ0P/`, `/PBLOCKMOMJ0M/`, `/PBLOCKMOMJ1/`, `/PBLOCKMOMJ2P/`, `/PBLOCKMOMJ2M/`, `/PBLOCKMOMJ0/`, `/PBLOCKMOMJ2/`, `/PBLOCKMOMJALL/`, `/LINES/`, `/LINESMOM/`, `/ITEM/`

### `SETUPB`

- Location: `torelons_4D_adjoint.f:17678`
- Kind/signature: `SUBROUTINE SETUPB`
- Description: Builds neighbor-index lookup tables for the blocked spatial lattice at each blocking level.
- Calls: none
- Called by: `MEASURE`
- Shared state: `/NEXTB/`

### `ACTION`

- Location: `torelons_4D_adjoint.f:17755`
- Kind/signature: `SUBROUTINE ACTION(IPR,ITER,TOTACT)`
- Description: Computes plaquette/action observables from gauge links and stores binned action/plaquette summaries.
- Calls: `VMX`, `HERM`, `TRVMX`, `JACKM`
- Called by: `MAIN`
- Shared state: `/ASTORE/`, `/COMMUNICATE/`, `/ARRAYS/`, `/NEXT/`

### `POLY`

- Location: `torelons_4D_adjoint.f:17880`
- Kind/signature: `SUBROUTINE POLY(IPR,ITER)`
- Description: Older or optional Polyakov-line measurement routine that computes temporal loops and jackknife averages.
- Calls: `VMX`, `JACKM`
- Called by: none
- Shared state: `/POLYT/`, `/ARRAYS/`, `/COMMUNICATE2/`, `/NEXT/`

### `RENORM`

- Location: `torelons_4D_adjoint.f:17987`
- Kind/signature: `SUBROUTINE RENORM`
- Description: Renormalizes all gauge links in `/ARRAYS/` by determinant normalization.
- Calls: `DETNANT`
- Called by: none
- Shared state: `/ARRAYS/`

### `DETNANT`

- Location: `torelons_4D_adjoint.f:18063`
- Kind/signature: `SUBROUTINE DETNANT(NUMOP,DET,AA)`
- Description: Copies a complex matrix into a work array and computes its determinant through LU decomposition.
- Calls: `DETMAT`
- Called by: `DODET`, `THERML1`, `RENORM`
- Shared state: none

### `DETMAT`

- Location: `torelons_4D_adjoint.f:18080`
- Kind/signature: `SUBROUTINE DETMAT(A,NP,DET)`
- Description: Computes a determinant from an LU factorization by multiplying the diagonal pivots and parity factor.
- Calls: `LUDCMP`
- Called by: `DETNANT`
- Shared state: none

### `LUDCMP`

- Location: `torelons_4D_adjoint.f:18098`
- Kind/signature: `SUBROUTINE LUDCMP(A,N,NP,INDX,D)`
- Description: Performs LU decomposition with partial pivoting, adapted from the classic Numerical Recipes-style routine.
- Calls: none
- Called by: `DETMAT`
- Shared state: none

### `SETUP`

- Location: `torelons_4D_adjoint.f:18168`
- Kind/signature: `SUBROUTINE SETUP`
- Description: Builds nearest-neighbor index tables for the full four-dimensional lattice.
- Calls: none
- Called by: `MAIN`
- Shared state: `/NEXT/`

### `JACKM`

- Location: `torelons_4D_adjoint.f:18224`
- Kind/signature: `SUBROUTINE JACKM(NUM,AV,AVER,ERR)`
- Description: Computes jackknife mean and error for a one-dimensional sample array.
- Calls: none
- Called by: `MLINE1`, `MLINE2`, `MLINE3`, `MLINE4`, `MLINE5`, `ACTION`, `POLY`
- Shared state: none

### `JACKMR`

- Location: `torelons_4D_adjoint.f:18250`
- Kind/signature: `SUBROUTINE JACKMR(NUM,NMAX,AVU,AVD,AVER,ERR)`
- Description: Computes jackknife estimates for ratios or paired numerator/denominator sample arrays.
- Calls: none
- Called by: none
- Shared state: none

### `JACK`

- Location: `torelons_4D_adjoint.f:18296`
- Kind/signature: `SUBROUTINE JACK(NBIN,ISUB,LXI,NUMBIN,LMAX,COR,VAC)`
- Description: Jackknife-analyzes correlators over bins/time separations and fills fit-summary arrays.
- Calls: `FITT`
- Called by: `MLINE1`, `MLINE2`, `MLINE3`, `MLINE4`, `MLINE5`
- Shared state: `/FIT/`

### `FITT`

- Location: `torelons_4D_adjoint.f:18385`
- Kind/signature: `SUBROUTINE FITT(LT,AW,BW,NTMAX)`
- Description: Estimates effective masses or fit-like transformed correlator quantities from time-slice data.
- Calls: none
- Called by: `JACK`
- Shared state: none

### `VMX`

- Location: `torelons_4D_adjoint.f:18431`
- Kind/signature: `SUBROUTINE VMX(NNN1,A,B,C,NNN2)`
- Description: Multiplies two `NCOL x NCOL` complex matrices stored in flattened vector form.
- Calls: none
- Called by: `POLYLOOP`, `BLOCK`, `SMEAR1`, `DIAG`, `RENORMBS`, `THERML1`, `ACTION`, `POLY`
- Shared state: none

### `TRVMX`

- Location: `torelons_4D_adjoint.f:18454`
- Kind/signature: `SUBROUTINE TRVMX(NNN1,A,B,CC,NNN2)`
- Description: Multiplies two flattened complex matrices and returns the trace of the product.
- Calls: none
- Called by: `POLYLOOP`, `ACTION`
- Shared state: none

### `HERM`

- Location: `torelons_4D_adjoint.f:18474`
- Kind/signature: `SUBROUTINE HERM(NNN1,A11,DUM11,NNN2)`
- Description: Forms the Hermitian conjugate of a flattened complex matrix.
- Calls: none
- Called by: `SMEAR1`, `DIAG`, `THERML1`, `ACTION`
- Shared state: none

### `RNDNUM`

- Location: `luesher.f:47`
- Kind/signature: `real*8 function rndnum()`
- Description: Returns one double-precision random number using the Luescher/Wolff RCARRY generator state.
- Calls: none
- Called by: none
- Shared state: `/RRAND/`

### `MULTI_RNDNUM`

- Location: `luesher.f:107`
- Kind/signature: `SUBROUTINE MULTI_RNDNUM(num,rndnum)`
- Description: Fills an array with multiple random numbers from the same RCARRY generator to reduce call overhead.
- Calls: none
- Called by: none
- Shared state: `/RRAND/`

### `RINI`

- Location: `luesher.f:175`
- Kind/signature: `subroutine rini(iran)`
- Description: Initializes the RCARRY random number generator state from an integer seed.
- Calls: none
- Called by: `MAIN`
- Shared state: `/RRAND/`

## Shared `COMMON` Blocks

- `/ARRAYB/`: `BLOCK`, `THERML1`
- `/ARRAYS/`: `MAIN`, `READ_GF`, `READ_GAUGE_ETMC`, `POLYLOOP`, `BLOCK`, `ACTION`, `POLY`, `RENORM`
- `/ASMEAR1/`: `BLOCK`, `SMEAR1`, `DODET`
- `/ASMEAR2/`: `BLOCK`, `SMEAR1`
- `/ASTORE/`: `TODISK1`, `ACTION`
- `/COMMUNICATE/`: `ACTION`
- `/COMMUNICATE2/`: `POLY`
- `/DIAGIN/`: `SMEAR1`, `DIAG`
- `/DIAGOUT/`: `SMEAR1`, `DIAG`
- `/DUMMY/`: `THERML1`
- `/FIT/`: `MLINE1`, `MLINE2`, `MLINE3`, `MLINE4`, `MLINE5`, `JACK`
- `/ITEM/`: `MEASURE`, `MLINE1`, `MLINE2`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/LINES/`: `THERML1`, `POT`
- `/LINESMOM/`: `THERML1`, `POT`
- `/NEXT/`: `POLYLOOP`, `ACTION`, `POLY`, `SETUP`
- `/NEXTB/`: `BLOCK`, `THERML1`, `SETUPB`
- `/PBLOCKJ0/`: `TODISK1`, `POT`
- `/PBLOCKJ0MM/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ0MP/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ0PM/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ0PP/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ1/`: `TODISK1`, `POT`
- `/PBLOCKJ1M/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ1P/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ2/`: `TODISK1`, `POT`
- `/PBLOCKJ2MM/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ2MP/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ2PM/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJ2PP/`: `TODISK1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKJALL/`: `TODISK1`, `POT`
- `/PBLOCKMOMJ0/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ0M/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ0P/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ1/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ2/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ2M/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJ2P/`: `TODISK2`, `MLINE4`, `MLINE5`, `POT`
- `/PBLOCKMOMJALL/`: `TODISK2`, `POT`
- `/PBLOCKMOMP/`: `TODISK2`, `MLINE2`, `POT`
- `/PBLOCKP/`: `TODISK1`, `MLINE1`, `MLINE3`, `MLINE4`, `MLINE5`, `POT`
- `/POLYT/`: `TODISK1`, `POLY`
- `/RRAND/`: `RNDNUM`, `MULTI_RNDNUM`, `RINI`
- `/SUBB/`: `RUNGRB`, `SUBGRB`

## Notes

- Defined but not reached from the current `MAIN` call tree: `QUATERNION_TO_MATRIX`, `READ_GAUGE_ETMC`, `POLY`, `RENORM`, `JACKMR`, `RNDNUM`, `MULTI_RNDNUM`.
- `READ_GF` declares `quaternion_to_matrix`, but the routine currently fills the matrix inline instead of calling it.
- Function-style dependencies are intentionally not inferred from `name(...)` text because this fixed-form code also uses array syntax with the same shape; explicit `CALL` edges are listed exactly.
- This is a static scan, so dependencies hidden behind generated code, compiler-specific behavior, or dynamic command strings are not expanded beyond the visible source.
