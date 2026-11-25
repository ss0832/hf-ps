# hf-ps

A minimal, pure PowerShell implementation of the restricted Hartree–Fock (RHF) method. It computes molecular integrals and solves the self-consistent field (SCF) equations without external numerical libraries, compiled extensions, or Python at runtime.

**Version:** 0.1.0 (experimental)

The initial implementation targets the hydrogen molecule (H2) in the STO-3G basis. It is a small electronic-structure codebase rather than a general-purpose quantum chemistry package.

## Features

- Contracted, normalized primitive s-type Gaussian basis functions
- Analytical overlap, kinetic-energy, nuclear-attraction, and electron-repulsion integrals
- Numerically evaluated zeroth-order Boys function
- Closed-shell Fock matrix and density matrix construction
- Symmetry-adapted solution of the Roothaan–Hall equations for H2
- Energy and density convergence criteria
- Molecular orbital energies, SCF history, and reusable integral tensors
- Command-line interface, PowerShell module API, bond scans, and regression tests

## Requirements

- Windows PowerShell 5.1 or PowerShell 7 (target compatibility)
- No additional packages or executables

The original standalone script was successfully executed on Windows PowerShell at a bond length of 1.4 bohr. The modular refactor has been statically reviewed, but it has not yet been executed in a PowerShell runtime. The included CI workflow and test script are intended to validate this refactor.

## Quick start

Open PowerShell in the repository root and run:

```powershell
.\hf-ps.ps1
```

Set the H–H distance in bohr and the SCF settings:

```powershell
.\hf-ps.ps1 -BondLength 1.4 -MaxIter 100 -Tolerance 1e-10
```

If local execution policy blocks a downloaded script, inspect its contents first. On a system you administer, a session-only execution policy can be used if appropriate:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
.\hf-ps.ps1
```

This setting does not persist after the session ends, and it does not override organization-level policy.

## Reference output

The following output was obtained with the original standalone implementation included under `examples/`. It is the expected regression target for the modular CLI.

```text
H2 RHF/STO-3G  R = 1.400000 bohr
 Iter        E(total) / Eh           |dE|        max|dP|
    1      -1.116714325063              ∞    6.0266E-001
    2      -1.116714325063    0.0000E+000    0.0000E+000

SCF converged.
Electronic energy  = -1.831000039348 Eh
Nuclear repulsion  = 0.714285714286 Eh
Total energy       = -1.116714325063 Eh
sigma_g energy     = -0.578202977512 Eh
sigma_u energy     = 0.670267768274 Eh
Tr(P S)            = 2.000000000000
S(1,2)             = 0.659318200128
```

Formatting of infinity and decimal numbers may vary with the PowerShell version and system locale.

## Module API

```powershell
Import-Module ./src/hf-ps.psd1 -Force

$result = Invoke-HfRHF -BondLength 1.4
$result.TotalEnergy
$result.SigmaGEnergy
$result.SigmaUEnergy
$result.Density[0,1]
$result.History | Format-Table

$integrals = Get-HfH2Integrals -BondLength 1.4
$integrals.Overlap[0,1]
$integrals.ElectronRepulsion[0,0,1,1]
```

The command-line script also supports structured output:

```powershell
$result = .\hf-ps.ps1 -AsObject
$result.TotalEnergy
```

All quantities use atomic units: coordinates and bond lengths in bohr; energies in hartree. PowerShell array indices are zero-based. Electron-repulsion integrals use chemists' notation `(ij|kl)`.

### Public functions

| Function | Purpose |
| --- | --- |
| `Get-HfH2Integrals` | Calculate the AO overlap, kinetic, nuclear-attraction, core Hamiltonian, and four-index ERI tensors. |
| `Invoke-HfRHF` | Run the H2 RHF calculation and return a structured result containing energies, orbital coefficients, density, Fock matrix, integrals, and SCF history. |

Both functions accept `-BondLength` (default 1.4 bohr). `Invoke-HfRHF` also accepts `-MaxIter` (default 100) and `-Tolerance` (default 1e-10), used for both absolute energy change and maximum absolute density-matrix change. It throws an error if SCF does not converge.

### Bond-length scan

```powershell
.\examples\bond-scan.ps1 | Format-Table
.\examples\bond-scan.ps1 -OutputPath ./h2-scan.csv
```

The scan uses the same RHF solver and does not perform geometry optimization or correlated calculations.

## Theory

The basis consists of two contracted 1s Gaussian functions, each containing three normalized primitives. The one-electron Hamiltonian is `H = T + V`. With a spin-summed closed-shell density `P`, the Fock matrix is

```text
F_ij = H_ij + sum_kl P_kl [(ij|kl) - 0.5 (ik|jl)]
P_ij = 2 C_i C_j
```

The generalized eigenvalue problem is `F C = S C epsilon`. For this symmetric two-center, two-AO problem, the gerade and ungerade vectors are proportional to `(1, 1)` and `(1, -1)`, respectively, normalized using the overlap matrix. This avoids a general eigensolver without replacing the integral evaluation or Fock construction.

The electronic and total energies are

```text
E_elec  = 0.5 sum_ij P_ij (H_ij + F_ij)
E_total = E_elec + 1/R
```

The energy is evaluated using a Fock matrix built from the updated density. The default H2 calculation converges in two iterations because the symmetric occupied orbital is fixed by the two-center symmetry.

## Project layout

```text
hf-ps/
  hf-ps.ps1                        CLI entry point
  src/
    hf-ps.psd1                     Module manifest
    hf-ps.psm1                     Module loader
    Private/
      Basis.ps1                    Hydrogen STO-3G parameters
      Boys.ps1                     F0 Boys function
      Fock.ps1                     RHF Fock builder
    Public/
      Get-HfH2Integrals.ps1         Analytical integral engine
      Invoke-HfRHF.ps1             SCF solver and result object
  examples/
    bond-scan.ps1                  Bond-length scan
    h2-rhf-standalone.ps1           Original user-validated script
  tests/
    run-tests.ps1                  Dependency-free regression tests
  .github/workflows/tests.yml       CI configuration
  README.md
  CHANGELOG.md
  .gitignore
```

## Testing

```powershell
.\tests\run-tests.ps1
```

No Pester installation is required. The tests cover the original Windows PowerShell reference values, several other bond lengths, AO/ERI symmetry, density idempotency, the Roothaan–Hall residual, orbital normalization, energy decomposition, public API and CLI outputs, and nonconvergence handling. The additional bond-length energy targets were calculated independently with an error-function-based Boys evaluator.

The GitHub Actions workflow defines tests on Windows PowerShell 5.1, PowerShell 7 on Windows, and PowerShell 7 on Linux. These jobs have not yet been run.

## Scope and limitations

- H2 only, with two 1s STO-3G basis functions and two electrons.
- Restricted, closed-shell Hartree–Fock; no UHF, ROHF, DFT, or post-HF methods.
- No general molecular input, angular-momentum recurrences, basis-set parser, or generalized eigensolver.
- No analytical gradients, geometry optimization, or stability analysis.
- Restricted HF does not correctly describe the H2 dissociation limit. Near-degenerate orbitals at very large distances can also affect numerical behavior. The accepted input range is not a validated accuracy range.
- The implementation prioritizes transparency over computational performance.

## Possible extensions

General molecular geometry and basis representation, higher-angular-momentum integrals, symmetric orthogonalization and a Jacobi eigensolver, arbitrary closed-shell molecules, and numerical or analytical energy derivatives can be developed independently of the CLI.

## Development note

This codebase was generated with AI assistance. The standalone baseline has been user-executed; the modular implementation still requires runtime testing. No license has been selected for this repository.
