param(
    [ValidateRange(0.05, 100.0)]
    [double]$BondLength = 1.4,
    [ValidateRange(1, 10000)]
    [int]$MaxIter = 100,
    [ValidateRange(1e-15, 1.0)]
    [double]$Tolerance = 1e-10
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# H2, atomic units, normalized primitive 1s STO-3G basis.
$pi = [Math]::PI
[double[]]$alpha = @(3.42525091, 0.62391373, 0.16885540)
[double[]]$coeff = @(0.15432897, 0.53532814, 0.44463454)
[double[]]$centers = @(0.0, $BondLength)

function Get-Boys0([double]$t) {
    # F0(t) = integral_0^1 exp(-t*u*u) du.
    if ($t -gt 40.0) {
        return (0.5 * [Math]::Sqrt([Math]::PI / $t))
    }
    # Positive series: exp(-t) * sum_k (2t)^k/(2k+1)!!
    $term = 1.0
    $sum = 1.0
    for ($n = 1; $n -lt 512; $n++) {
        $term *= 2.0 * $t / (2.0 * $n + 1.0)
        $sum += $term
        if ($term -le 1e-16 * $sum) { break }
    }
    return ([Math]::Exp(-$t) * $sum)
}

# Precompute Gaussian product data for each AO pair.
$weights = [double[]]::new(3)
for ($a = 0; $a -lt 3; $a++) {
    $weights[$a] = $coeff[$a] * [Math]::Pow(2.0 * $alpha[$a] / $pi, 0.75)
}
$pairs = [object[,]]::new(2, 2)
for ($i = 0; $i -lt 2; $i++) {
    for ($j = 0; $j -lt 2; $j++) {
        $dz = $centers[$i] - $centers[$j]
        $r2 = $dz * $dz
        $pairs[$i,$j] = @(
            for ($a = 0; $a -lt 3; $a++) {
                for ($b = 0; $b -lt 3; $b++) {
                    $p = $alpha[$a] + $alpha[$b]
                    $mu = $alpha[$a] * $alpha[$b] / $p
                    [pscustomobject]@{
                        Exponent = $p
                        Z = ($alpha[$a] * $centers[$i] + $alpha[$b] * $centers[$j]) / $p
                        MuR2 = $mu * $r2
                        Mu = $mu
                        Weight = $weights[$a] * $weights[$b] * [Math]::Exp(-$mu * $r2)
                    }
                }
            }
        )
    }
}

# One-electron integrals: overlap S, kinetic T, attraction V, core H.
$S = [double[,]]::new(2, 2)
$T = [double[,]]::new(2, 2)
$V = [double[,]]::new(2, 2)
$H = [double[,]]::new(2, 2)
for ($i = 0; $i -lt 2; $i++) {
    for ($j = 0; $j -lt 2; $j++) {
        foreach ($ab in $pairs[$i,$j]) {
            $p = $ab.Exponent
            $ov = $ab.Weight * [Math]::Pow($pi / $p, 1.5)
            $S[$i,$j] += $ov
            $T[$i,$j] += $ab.Mu * (3.0 - 2.0 * $ab.MuR2) * $ov
            foreach ($zc in $centers) {
                $d = $ab.Z - $zc
                $f0 = Get-Boys0 ($p * $d * $d)
                $V[$i,$j] -= $ab.Weight * (2.0 * $pi / $p) * $f0
            }
        }
        $H[$i,$j] = ($T[$i,$j]) + ($V[$i,$j])
    }
}

# Electron repulsion integrals, chemists' notation (ij|kl).
$ERI = [double[,,,]]::new(2, 2, 2, 2)
$eriPref = 2.0 * [Math]::Pow($pi, 2.5)
for ($i = 0; $i -lt 2; $i++) {
    for ($j = 0; $j -lt 2; $j++) {
        for ($k = 0; $k -lt 2; $k++) {
            for ($l = 0; $l -lt 2; $l++) {
                foreach ($ab in $pairs[$i,$j]) {
                    foreach ($cd in $pairs[$k,$l]) {
                        $p = $ab.Exponent
                        $q = $cd.Exponent
                        $d = $ab.Z - $cd.Z
                        $f0 = Get-Boys0 ($p * $q / ($p + $q) * $d * $d)
                        $ERI[$i,$j,$k,$l] += $ab.Weight * $cd.Weight * $eriPref * $f0 / ($p * $q * [Math]::Sqrt($p + $q))
                    }
                }
            }
        }
    }
}

function Get-Fock([double[,]]$density) {
    $fock = [double[,]]::new(2, 2)
    for ($i = 0; $i -lt 2; $i++) {
        for ($j = 0; $j -lt 2; $j++) {
            $fock[$i,$j] = $H[$i,$j]
            for ($k = 0; $k -lt 2; $k++) {
                for ($l = 0; $l -lt 2; $l++) {
                    $coulomb = $ERI[$i,$j,$k,$l]
                    $exchange = $ERI[$i,$k,$j,$l]
                    $fock[$i,$j] += ($density[$k,$l]) * ($coulomb - 0.5 * $exchange)
                }
            }
        }
    }
    # Preserve the multidimensional array as a single pipeline object.
    return ,$fock
}

# SCF: H2 symmetry gives the S-normalized gerade/ungerade MOs.
$D = [double[,]]::new(2, 2)
$sPlus = 0.5 * (($S[0,0]) + ($S[1,1]) + ($S[0,1]) + ($S[1,0]))
$sMinus = 0.5 * (($S[0,0]) + ($S[1,1]) - ($S[0,1]) - ($S[1,0]))
$eNuc = 1.0 / $BondLength
$ePrev = [double]::PositiveInfinity
$converged = $false
Write-Host ('H2 RHF/STO-3G  R = {0:F6} bohr' -f $BondLength)
Write-Host ('{0,5} {1,20} {2,14} {3,14}' -f 'Iter','E(total) / Eh','|dE|','max|dP|')
for ($iter = 1; $iter -le $MaxIter; $iter++) {
    $F = Get-Fock $D
    $fDiag = 0.5 * (($F[0,0]) + ($F[1,1]))
    $fOff = 0.5 * (($F[0,1]) + ($F[1,0]))
    $epsG = ($fDiag + $fOff) / $sPlus
    $epsU = ($fDiag - $fOff) / $sMinus
    if ($epsG -le $epsU) {
        $c0 = 1.0 / [Math]::Sqrt(2.0 * $sPlus)
        $c1 = $c0
    } else {
        $c0 = 1.0 / [Math]::Sqrt(2.0 * $sMinus)
        $c1 = -$c0
    }
    [double[]]$C = @($c0, $c1)
    $Dnew = [double[,]]::new(2, 2)
    $dP = 0.0
    for ($i = 0; $i -lt 2; $i++) {
        for ($j = 0; $j -lt 2; $j++) {
            $Dnew[$i,$j] = 2.0 * ($C[$i]) * ($C[$j])
            $delta = ($Dnew[$i,$j]) - ($D[$i,$j])
            $dP = [Math]::Max($dP, [Math]::Abs($delta))
        }
    }
    # Use F(Dnew) for a consistent energy at every iteration.
    $Fnew = Get-Fock $Dnew
    $eElec = 0.0
    for ($i = 0; $i -lt 2; $i++) {
        for ($j = 0; $j -lt 2; $j++) {
            $eElec += 0.5 * ($Dnew[$i,$j]) * (($H[$i,$j]) + ($Fnew[$i,$j]))
        }
    }
    $eTotal = $eElec + $eNuc
    $dE = [Math]::Abs($eTotal - $ePrev)
    Write-Host ('{0,5} {1,20:F12} {2,14:E4} {3,14:E4}' -f $iter,$eTotal,$dE,$dP)
    $D = $Dnew
    if (($dE -lt $Tolerance) -and ($dP -lt $Tolerance)) {
        $converged = $true
        break
    }
    $ePrev = $eTotal
}
if (-not $converged) { throw 'SCF did not converge.' }

# Final eigenvalues, evaluated from the converged density.
$fDiag = 0.5 * (($Fnew[0,0]) + ($Fnew[1,1]))
$fOff = 0.5 * (($Fnew[0,1]) + ($Fnew[1,0]))
$epsG = ($fDiag + $fOff) / $sPlus
$epsU = ($fDiag - $fOff) / $sMinus
$electronCount = 0.0
for ($i = 0; $i -lt 2; $i++) {
    for ($j = 0; $j -lt 2; $j++) {
        $electronCount += ($D[$i,$j]) * ($S[$j,$i])
    }
}
Write-Host ''
Write-Host 'SCF converged.'
Write-Host ('Electronic energy  = {0:F12} Eh' -f $eElec)
Write-Host ('Nuclear repulsion  = {0:F12} Eh' -f $eNuc)
Write-Host ('Total energy       = {0:F12} Eh' -f $eTotal)
Write-Host ('sigma_g energy     = {0:F12} Eh' -f $epsG)
Write-Host ('sigma_u energy     = {0:F12} Eh' -f $epsU)
Write-Host ('Tr(P S)            = {0:F12}' -f $electronCount)
Write-Host ('S(1,2)             = {0:F12}' -f ($S[0,1]))
