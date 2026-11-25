function Get-HfHydrogenSto3g {
    # Normalized primitive 1s Gaussian functions, atomic units.
    [double[]]$alpha = @(3.42525091, 0.62391373, 0.16885540)
    [double[]]$coeff = @(0.15432897, 0.53532814, 0.44463454)
    $weights = [double[]]::new(3)
    for ($a = 0; $a -lt 3; $a++) {
        $weights[$a] = $coeff[$a] * [Math]::Pow(2.0 * $alpha[$a] / [Math]::PI, 0.75)
    }
    return [pscustomobject]@{
        Exponents = $alpha
        ContractionCoefficients = $coeff
        Weights = $weights
    }
}
