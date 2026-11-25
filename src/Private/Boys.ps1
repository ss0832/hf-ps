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
