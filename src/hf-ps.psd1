@{
    RootModule = 'hf-ps.psm1'
    ModuleVersion = '0.1.0'
    GUID = '7a2464aa-5e25-428d-9906-e836e6f24711'
    Description = 'A pure PowerShell restricted Hartree-Fock implementation for H2.'
    PowerShellVersion = '5.1'
    FunctionsToExport = @('Get-HfH2Integrals', 'Invoke-HfRHF')
    CmdletsToExport = @()
    VariablesToExport = @()
    AliasesToExport = @()
}
