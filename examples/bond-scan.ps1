[CmdletBinding()]
param(
    [double[]]$BondLengths = @(0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0),
    [string]$OutputPath
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Import-Module (Join-Path (Split-Path -Parent $PSScriptRoot) 'src/hf-ps.psd1') -Force
$rows = @(
    foreach ($r in $BondLengths) {
        $result = Invoke-HfRHF -BondLength $r
        [pscustomobject]@{
            BondLengthBohr = $r
            TotalEnergyEh = $result.TotalEnergy
            ElectronicEnergyEh = $result.ElectronicEnergy
            NuclearRepulsionEh = $result.NuclearRepulsion
        }
    }
)
if ($OutputPath) {
    $rows | Export-Csv -Path $OutputPath -NoTypeInformation -Encoding UTF8
} else {
    $rows
}
