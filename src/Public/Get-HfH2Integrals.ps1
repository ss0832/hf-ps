function Get-HfH2Integrals {
    [CmdletBinding()]
    param(
        [ValidateRange(0.05, 100.0)]
        [double]$BondLength = 1.4
    )

    $pi = [Math]::PI
    $basis = Get-HfHydrogenSto3g
    $alpha = $basis.Exponents
    $weights = $basis.Weights
    [double[]]$centers = @(0.0, $BondLength)

    # Precompute Gaussian product data for each AO pair.
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

    return [pscustomobject]@{
        BondLength = $BondLength
        Basis = 'STO-3G'
        Overlap = $S
        Kinetic = $T
        NuclearAttraction = $V
        CoreHamiltonian = $H
        ElectronRepulsion = $ERI
        NuclearRepulsion = 1.0 / $BondLength
    }
}
