# Builds and checks the package in a scratch copy, so that compiled objects and check
# output never land in the source folder (which sits in a synced directory).
#   powershell -File tools/build_check.ps1 <scratch_dir> [docs|install|check]
# On Linux/macOS: R CMD build . && R CMD check --as-cran PrecFact_*.tar.gz
param([Parameter(Mandatory = $true)][string]$B, [string]$Mode = "check")
$PKG = Split-Path -Parent $PSScriptRoot
if (Test-Path "$B\PrecFact") { Remove-Item -Recurse -Force "$B\PrecFact" }
New-Item -ItemType Directory -Force "$B\lib" | Out-Null
robocopy $PKG "$B\PrecFact" /E /XD .git figures /NFL /NDL /NJH /NJS | Out-Null
Set-Location $B
switch ($Mode) {
  "docs" {
    Set-Content -Path docs.R -Encoding ascii -Value 'Rcpp::compileAttributes("PrecFact"); roxygen2::roxygenise("PrecFact")'
    & Rscript docs.R
    Copy-Item "PrecFact\NAMESPACE" "$PKG\NAMESPACE" -Force
    Copy-Item "PrecFact\R\RcppExports.R" "$PKG\R\" -Force
    Copy-Item "PrecFact\src\RcppExports.cpp" "$PKG\src\" -Force
    Remove-Item "$PKG\man\*.Rd"; Copy-Item "PrecFact\man\*.Rd" "$PKG\man\" -Force
  }
  "install" { & R.exe CMD INSTALL --library="$B\lib" PrecFact }
  "check" {
    & R.exe CMD build PrecFact
    & R.exe CMD check --as-cran --no-manual (Get-ChildItem PrecFact_*.tar.gz).Name
  }
}
