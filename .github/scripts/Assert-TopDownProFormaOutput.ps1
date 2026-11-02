<#
.SYNOPSIS
    Asserts that a top-down search actually produced proteoforms and wrote a populated ProForma column.

.DESCRIPTION
    The ProForma column is emitted only for top-down runs (it is gated on AnalyteType.Proteoform),
    and it is the only thing in MetaMorpheus that calls into mzLib's Readers.ProForma. A search that
    finds zero proteoforms therefore never reaches ToProFormaString(), and will report success while
    leaving that whole path untested - which is how a TopDownProteomics version mismatch reached
    master green and only failed on a real search.

    Checking the exit code alone is not enough. This asserts the column exists AND carries a value,
    so the writer is proven to have run.
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string] $SearchOutput
)

$ErrorActionPreference = 'Stop'

$results = Get-ChildItem -Path $SearchOutput -Filter 'AllProteoforms.psmtsv' -Recurse -ErrorAction SilentlyContinue
if (-not $results) {
    Write-Host "AllProteoforms.psmtsv was not written under $SearchOutput" -ForegroundColor Red
    exit 1
}

$resultFile = $results[0].FullName
$lines = Get-Content -LiteralPath $resultFile
if ($lines.Count -lt 2) {
    Write-Host "$resultFile has no proteoform rows, so the ProForma writer never ran." -ForegroundColor Red
    Write-Host "This fixture is only meaningful if the search finds something - fix the fixture, not this check." -ForegroundColor Red
    exit 1
}

$header = $lines[0] -split "`t"
$proFormaIndex = [Array]::IndexOf($header, 'ProForma')
if ($proFormaIndex -lt 0) {
    Write-Host "No ProForma column in the header of $resultFile" -ForegroundColor Red
    exit 1
}

$firstRow = $lines[1] -split "`t"
if ($firstRow.Count -ne $header.Count) {
    Write-Host "Row/header column count mismatch ($($firstRow.Count) vs $($header.Count)) in $resultFile" -ForegroundColor Red
    exit 1
}

$proFormaValue = $firstRow[$proFormaIndex]
if ([string]::IsNullOrWhiteSpace($proFormaValue)) {
    Write-Host "ProForma column is empty in $resultFile" -ForegroundColor Red
    exit 1
}

Write-Host "Top-down search wrote $($lines.Count - 1) proteoform row(s) with a populated ProForma column." -ForegroundColor Green
Write-Host "  ProForma: $proFormaValue"
