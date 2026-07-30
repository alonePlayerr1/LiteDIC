<#
  sync_engine_to_app.ps1
  Keep the LiteDIC app's vendored engine\ folder (6_DIC + 11_Stereo, *.m only)
  and this lab checkout in agreement.

  Three modes:

      -Check          report drift both ways, change nothing.
                      Exit code 1 when anything differs, so it can gate a
                      commit or a test run.
      (default)       PULL  lab -> app. Refuses to overwrite app-side files
                      that differ from the lab, because those are patches
                      nobody has mirrored back yet; -Force overrides.
      -Push           PUSH  app -> lab. The direction that keeps happening in
                      practice: a defect is found while running the app, fixed
                      in engine\, and then copied back by hand.

  Run from anywhere:
      powershell -ExecutionPolicy Bypass -File "sync_engine_to_app.ps1" -Check
  Override the app location if it is not the sibling "LiteDIC app":
      ... -File sync_engine_to_app.ps1 -App "D:\path\to\LiteDIC app"

  ENGINE_SOURCE.txt is re-stamped, not rewritten: the 'patched:' notes are the
  only record of where and why the vendored copy diverged.
#>
param(
    [string]$App = (Join-Path (Split-Path $PSScriptRoot -Parent) 'LiteDIC app'),
    [switch]$Check,
    [switch]$Push,
    [switch]$Force
)
$ErrorActionPreference = 'Stop'
$lab = $PSScriptRoot
$dst = Join-Path $App 'engine'
$dirs = @('6_DIC', '11_Stereo')

if (-not (Test-Path $App)) { throw "App folder not found: $App" }
if ($Check -and $Push)     { throw "-Check and -Push are mutually exclusive." }

function Get-Drift {
    # Per module: which .m differ, which exist on one side only.
    $out = @()
    foreach ($f in $dirs) {
        $labDir = Join-Path $lab $f
        $appDir = Join-Path $dst $f
        $labFiles = @{}
        $appFiles = @{}
        if (Test-Path $labDir) {
            Get-ChildItem $labDir -Filter *.m | ForEach-Object { $labFiles[$_.Name] = $_.FullName }
        }
        if (Test-Path $appDir) {
            Get-ChildItem $appDir -Filter *.m | ForEach-Object { $appFiles[$_.Name] = $_.FullName }
        }
        foreach ($n in ($labFiles.Keys + $appFiles.Keys | Sort-Object -Unique)) {
            $inLab = $labFiles.ContainsKey($n)
            $inApp = $appFiles.ContainsKey($n)
            $state = ''
            if     (-not $inApp) { $state = 'lab only' }
            elseif (-not $inLab) { $state = 'app only' }
            else {
                $a = (Get-FileHash $labFiles[$n] -Algorithm SHA256).Hash
                $b = (Get-FileHash $appFiles[$n] -Algorithm SHA256).Hash
                if ($a -ne $b) { $state = 'DIFFERS' }
            }
            if ($state -ne '') {
                $out += [pscustomobject]@{ Module = $f; File = $n; State = $state }
            }
        }
    }
    return $out
}

function Write-Stamp {
    # Replace only the three provenance lines; keep every 'patched:' note and
    # anything else the file carries.
    $path = Join-Path $dst 'ENGINE_SOURCE.txt'
    $commit = (& git -C $lab rev-parse HEAD 2>$null)
    if (-not $commit) { $commit = '(git not available)' }
    $when = (Get-Date).ToUniversalTime().ToString('u') + ' UTC'

    if (Test-Path $path) {
        $lines = @(Get-Content $path)
        $kept = foreach ($l in $lines) {
            if ($l -match '^source_commit:') { "source_commit: $commit" }
            elseif ($l -match '^synced:')    { "synced:        $when" }
            else                             { $l }
        }
        if (-not ($kept -match '^source_commit:')) { $kept += "source_commit: $commit" }
        if (-not ($kept -match '^synced:'))        { $kept += "synced:        $when" }
        $nPatched = @($kept | Where-Object { $_ -match '^patched:' }).Count
        Set-Content -Path $path -Value $kept -Encoding UTF8
        Write-Host ("stamped ENGINE_SOURCE.txt (kept {0} 'patched:' note line(s))" -f $nPatched)
    } else {
        $stamp = @(
            'Vendored copy of the LiteDIC engine (code only: 6_DIC + 11_Stereo *.m).',
            '',
            'source_repo:   LiteDIC (lab)',
            "source_commit: $commit",
            "synced:        $when"
        )
        Set-Content -Path $path -Value $stamp -Encoding UTF8
        Write-Host 'wrote a new ENGINE_SOURCE.txt'
    }
}

# ---------------------------------------------------------------- check ----
$drift = @(Get-Drift)
if ($Check) {
    Write-Host "lab: $lab"
    Write-Host "app: $dst`n"
    if ($drift.Count -eq 0) {
        Write-Host 'IN SYNC - every vendored .m is byte-identical to the lab.'
        exit 0
    }
    $drift | Format-Table -AutoSize | Out-String | Write-Host
    Write-Host ("{0} file(s) differ. Use -Push (app -> lab) or the default pull (lab -> app)." -f $drift.Count)
    exit 1
}

# ----------------------------------------------------------------- push ----
if ($Push) {
    foreach ($f in $dirs) {
        $srcDir = Join-Path $dst $f
        $dstDir = Join-Path $lab $f
        if (-not (Test-Path $srcDir)) { continue }
        New-Item -ItemType Directory -Force -Path $dstDir | Out-Null
        Copy-Item (Join-Path $srcDir '*.m') $dstDir -Force
        $n = (Get-ChildItem $srcDir -Filter *.m).Count
        Write-Host ("pushed {0,-10} {1} files -> lab" -f $f, $n)
    }
    $left = @(Get-Drift)
    if ($left.Count -gt 0) {
        $left | Format-Table -AutoSize | Out-String | Write-Host
        throw "push finished but $($left.Count) file(s) still differ."
    }
    Write-Host "`nverified: every vendored .m is byte-identical to the lab."
    Write-Host 'NOT committed in the lab - run RUN_TESTS.m there, commit, then'
    Write-Host 're-run this script without -Push to re-stamp source_commit.'
    exit 0
}

# ----------------------------------------------------------------- pull ----
# App-side files that differ are unmirrored patches. Don't overwrite them
# silently, stop instead.
$risky = @($drift | Where-Object { $_.State -eq 'DIFFERS' -or $_.State -eq 'app only' })
if ($risky.Count -gt 0 -and -not $Force) {
    $risky | Format-Table -AutoSize | Out-String | Write-Host
    Write-Host 'These app-side files are NOT in the lab as they stand.'
    Write-Host 'Mirror them first:   sync_engine_to_app.ps1 -Push'
    Write-Host 'Or discard them:     sync_engine_to_app.ps1 -Force'
    throw "refusing to overwrite $($risky.Count) unmirrored app-side file(s)."
}

foreach ($f in $dirs) {
    $srcDir = Join-Path $lab $f
    $dstDir = Join-Path $dst $f
    New-Item -ItemType Directory -Force -Path $dstDir | Out-Null
    # drop stale .m that no longer exist in the lab, then mirror the code
    Get-ChildItem $dstDir -Filter *.m -ErrorAction SilentlyContinue |
        Where-Object { -not (Test-Path (Join-Path $srcDir $_.Name)) } |
        Remove-Item -Force
    Copy-Item (Join-Path $srcDir '*.m') $dstDir -Force
    $n = (Get-ChildItem $dstDir -Filter *.m).Count
    Write-Host ("synced {0,-10} {1} files" -f $f, $n)
}

Write-Stamp
$left = @(Get-Drift)
if ($left.Count -gt 0) {
    $left | Format-Table -AutoSize | Out-String | Write-Host
    throw "sync finished but $($left.Count) file(s) still differ."
}
Write-Host "verified in sync -> $dst"
