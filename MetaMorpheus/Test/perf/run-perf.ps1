# Runs the custom-oxonium performance bounds. These tests are [Explicit] + [Category("Performance")],
# so they are skipped by the normal `dotnet test` run and only execute when selected here.
#
# Usage:  pwsh Test\perf\run-perf.ps1     (from the MetaMorpheus solution directory)
$ErrorActionPreference = 'Stop'
$solutionDir = Resolve-Path (Join-Path $PSScriptRoot '..' '..')
Set-Location $solutionDir
dotnet test Test\Test.csproj --filter "Category=Performance" --nologo -v n
