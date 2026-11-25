[CmdletBinding()]
param()
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Import-Module (Join-Path $root 'src/hf-ps.psd1') -Force

function Assert-Near {
    param([string]$Name, [double]$Actual, [double]$Expected, [double]$Tolerance = 1e-10)
    if ([double]::IsNaN($Actual) -or [double]::IsInfinity($Actual) -or
        ([Math]::Abs($Actual - $Expected) -gt $Tolerance)) {
        throw ('{0}: expected {1:G17}, got {2:G17} (tolerance {3:G})' -f $Name,$Expected,$Actual,$Tolerance)
    }
}

# Regression against the standalone calculation executed on Windows PowerShell.
$r = Invoke-HfRHF -BondLength 1.4
if (-not $r.Converged) { throw 'SCF convergence failed.' }
Assert-Near 'SCF iterations' $r.Iterations 2 0
Assert-Near 'Total energy' $r.TotalEnergy -1.116714325063
Assert-Near 'Electronic energy' $r.ElectronicEnergy -1.831000039348
Assert-Near 'Nuclear repulsion' $r.NuclearRepulsion 0.714285714286
Assert-Near 'Sigma_g eigenvalue' $r.SigmaGEnergy -0.578202977512
Assert-Near 'Sigma_u eigenvalue' $r.SigmaUEnergy 0.670267768274
Assert-Near 'Electron count' $r.ElectronCount 2.0
Assert-Near 'AO overlap' $r.Overlap12 0.659318200128
Assert-Near 'Final energy change' $r.History[-1].EnergyChange 0.0
Assert-Near 'Final density change' $r.History[-1].DensityChange 0.0
Assert-Near 'Energy decomposition' $r.TotalEnergy ($r.ElectronicEnergy + $r.NuclearRepulsion)
Write-Host 'PASS: H2 reference calculation'

# Matrix symmetry, normalization, ERI permutation symmetry, and RHF identities.
$ints = $r.Integrals
$S = $ints.Overlap
$T = $ints.Kinetic
$V = $ints.NuclearAttraction
$H = $ints.CoreHamiltonian
$eri = $ints.ElectronRepulsion
$P = $r.Density
$F = $r.Fock
$C = $r.OccupiedCoefficients
Assert-Near 'AO normalization A' ($S[0,0]) 1.0 1e-7
Assert-Near 'AO normalization B' ($S[1,1]) 1.0 1e-7
$normalization = 0.0
$trace = 0.0
for ($i = 0; $i -lt 2; $i++) {
    $fc = 0.0
    $sc = 0.0
    for ($j = 0; $j -lt 2; $j++) {
        Assert-Near 'Overlap symmetry' ($S[$i,$j]) ($S[$j,$i])
        Assert-Near 'Kinetic symmetry' ($T[$i,$j]) ($T[$j,$i])
        Assert-Near 'Attraction symmetry' ($V[$i,$j]) ($V[$j,$i])
        Assert-Near 'Fock symmetry' ($F[$i,$j]) ($F[$j,$i])
        Assert-Near 'Core Hamiltonian' ($H[$i,$j]) (($T[$i,$j]) + ($V[$i,$j]))
        $normalization += ($C[$i]) * ($S[$i,$j]) * ($C[$j])
        $trace += ($P[$i,$j]) * ($S[$j,$i])
        $fc += ($F[$i,$j]) * ($C[$j])
        $sc += ($S[$i,$j]) * ($C[$j])
        $psp = 0.0
        for ($k = 0; $k -lt 2; $k++) {
            for ($l = 0; $l -lt 2; $l++) {
                $x = $eri[$i,$j,$k,$l]
                Assert-Near 'ERI first-pair permutation' $x ($eri[$j,$i,$k,$l])
                Assert-Near 'ERI second-pair permutation' $x ($eri[$i,$j,$l,$k])
                Assert-Near 'ERI pair exchange' $x ($eri[$k,$l,$i,$j])
                $psp += ($P[$i,$k]) * ($S[$k,$l]) * ($P[$l,$j])
            }
        }
        Assert-Near 'Closed-shell density idempotency' $psp (2.0 * ($P[$i,$j]))
    }
    Assert-Near 'Roothaan-Hall residual' $fc ($r.SigmaGEnergy * $sc)
}
Assert-Near 'Occupied MO normalization' $normalization 1.0
Assert-Near 'Independent electron count' $trace 2.0
Write-Host 'PASS: integral symmetry and RHF identities'

# Independent reference energies computed with an erf-based Boys function.
$cases = @(
    @{ R = 0.8; E = -0.947307931116 },
    @{ R = 1.0; E = -1.065999462143 },
    @{ R = 2.0; E = -1.049170901986 },
    @{ R = 4.0; E = -0.761082247024 }
)
foreach ($case in $cases) {
    $result = Invoke-HfRHF -BondLength $case.R
    Assert-Near ('Total energy at R={0} bohr' -f $case.R) $result.TotalEnergy $case.E
    Assert-Near 'Electron count along bond scan' $result.ElectronCount 2.0
}
Write-Host 'PASS: bond-length regression'

# Public integral API and object-mode CLI should give structured outputs.
$api = Get-HfH2Integrals
Assert-Near 'Integral API overlap' ($api.Overlap[0,1]) 0.659318200128
$cli = & (Join-Path $root 'hf-ps.ps1') -AsObject
Assert-Near 'CLI object-mode energy' $cli.TotalEnergy -1.116714325063
Write-Host 'PASS: public API and CLI'

$failedAsExpected = $false
try {
    $null = Invoke-HfRHF -MaxIter 1
} catch {
    if ($_.Exception.Message -notlike '*SCF did not converge*') { throw }
    $failedAsExpected = $true
}
if (-not $failedAsExpected) { throw 'Expected a nonconvergence exception.' }
Write-Host 'PASS: nonconvergence handling'
Write-Host 'All tests passed.'
