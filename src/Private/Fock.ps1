function Get-HfFock {
    param(
        [double[,]]$Density,
        [double[,]]$CoreHamiltonian,
        [double[,,,]]$ElectronRepulsion
    )
    $fock = [double[,]]::new(2, 2)
    for ($i = 0; $i -lt 2; $i++) {
        for ($j = 0; $j -lt 2; $j++) {
            $fock[$i,$j] = $CoreHamiltonian[$i,$j]
            for ($k = 0; $k -lt 2; $k++) {
                for ($l = 0; $l -lt 2; $l++) {
                    $coulomb = $ElectronRepulsion[$i,$j,$k,$l]
                    $exchange = $ElectronRepulsion[$i,$k,$j,$l]
                    $fock[$i,$j] += ($density[$k,$l]) * ($coulomb - 0.5 * $exchange)
                }
            }
        }
    }
    # Preserve the multidimensional array as a single pipeline object.
    return ,$fock
}
