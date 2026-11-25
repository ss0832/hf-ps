function Invoke-HfRHF {
    [CmdletBinding()]
    param(
        [ValidateRange(0.05, 100.0)]
        [double]$BondLength = 1.4,
        [ValidateRange(1, 10000)]
        [int]$MaxIter = 100,
        [ValidateRange(1e-15, 1.0)]
        [double]$Tolerance = 1e-10
    )

    $integrals = Get-HfH2Integrals -BondLength $BondLength
    $S = $integrals.Overlap
    $H = $integrals.CoreHamiltonian
    $ERI = $integrals.ElectronRepulsion

    # SCF: H2 symmetry gives the S-normalized gerade/ungerade MOs.
    $D = [double[,]]::new(2, 2)
    $sPlus = 0.5 * (($S[0,0]) + ($S[1,1]) + ($S[0,1]) + ($S[1,0]))
    $sMinus = 0.5 * (($S[0,0]) + ($S[1,1]) - ($S[0,1]) - ($S[1,0]))
    $eNuc = $integrals.NuclearRepulsion
    $ePrev = [double]::PositiveInfinity
    $converged = $false
    $history = @()
    for ($iter = 1; $iter -le $MaxIter; $iter++) {
        $F = Get-HfFock -Density $D -CoreHamiltonian $H -ElectronRepulsion $ERI
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
        $Fnew = Get-HfFock -Density $Dnew -CoreHamiltonian $H -ElectronRepulsion $ERI
        $eElec = 0.0
        for ($i = 0; $i -lt 2; $i++) {
            for ($j = 0; $j -lt 2; $j++) {
                $eElec += 0.5 * ($Dnew[$i,$j]) * (($H[$i,$j]) + ($Fnew[$i,$j]))
            }
        }
        $eTotal = $eElec + $eNuc
        $dE = [Math]::Abs($eTotal - $ePrev)
        $history += [pscustomobject]@{
            Iteration = $iter
            TotalEnergy = $eTotal
            EnergyChange = $dE
            DensityChange = $dP
        }
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

    return [pscustomobject]@{
        Molecule = 'H2'
        Method = 'RHF'
        Basis = 'STO-3G'
        BondLength = $BondLength
        Converged = $converged
        Iterations = $iter
        ElectronicEnergy = $eElec
        NuclearRepulsion = $eNuc
        TotalEnergy = $eTotal
        SigmaGEnergy = $epsG
        SigmaUEnergy = $epsU
        ElectronCount = $electronCount
        Overlap12 = $S[0,1]
        Density = $D
        Fock = $Fnew
        OccupiedCoefficients = $C
        History = $history
        Integrals = $integrals
    }
}
