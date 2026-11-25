[CmdletBinding()]
param(
    [ValidateRange(0.05, 100.0)]
    [double]$BondLength = 1.4,
    [ValidateRange(1, 10000)]
    [int]$MaxIter = 100,
    [ValidateRange(1e-15, 1.0)]
    [double]$Tolerance = 1e-10,
    [switch]$AsObject
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Import-Module (Join-Path $PSScriptRoot 'src/hf-ps.psd1') -Force
$result = Invoke-HfRHF -BondLength $BondLength -MaxIter $MaxIter -Tolerance $Tolerance
if ($AsObject) {
    return $result
}

Write-Host ('H2 RHF/STO-3G  R = {0:F6} bohr' -f $result.BondLength)
Write-Host ('{0,5} {1,20} {2,14} {3,14}' -f 'Iter','E(total) / Eh','|dE|','max|dP|')
foreach ($step in $result.History) {
    Write-Host ('{0,5} {1,20:F12} {2,14:E4} {3,14:E4}' -f $step.Iteration,$step.TotalEnergy,$step.EnergyChange,$step.DensityChange)
}
Write-Host ''
Write-Host 'SCF converged.'
Write-Host ('Electronic energy  = {0:F12} Eh' -f $result.ElectronicEnergy)
Write-Host ('Nuclear repulsion  = {0:F12} Eh' -f $result.NuclearRepulsion)
Write-Host ('Total energy       = {0:F12} Eh' -f $result.TotalEnergy)
Write-Host ('sigma_g energy     = {0:F12} Eh' -f $result.SigmaGEnergy)
Write-Host ('sigma_u energy     = {0:F12} Eh' -f $result.SigmaUEnergy)
Write-Host ('Tr(P S)            = {0:F12}' -f $result.ElectronCount)
Write-Host ('S(1,2)             = {0:F12}' -f $result.Overlap12)
