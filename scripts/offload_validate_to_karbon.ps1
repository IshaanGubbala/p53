param(
    [Parameter(Mandatory = $true)]
    [string]$RemoteUser,

    [Parameter(Mandatory = $true)]
    [string]$RemoteHost,

    [Parameter(Mandatory = $true)]
    [string]$RemoteProjectDir,

    [string]$RunId = "campaign_20260224_140210",

    [string]$RemotePython = "python",

    [string]$IdentityFile = "",

    [string]$ProxyJump = "",

    [string]$KexAlgorithms = "curve25519-sha256",

    [int]$ConnectTimeoutSec = 10,

    [bool]$BatchMode = $true,

    [switch]$IncludeEsmfold,

    [switch]$NoStopLocal,

    [switch]$AcceptNewHostKey
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Exe,
        [Parameter(Mandatory = $true)]
        [string[]]$Args
    )
    & $Exe @Args
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed ($LASTEXITCODE): $Exe $($Args -join ' ')"
    }
}

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$localCampaignDir = Join-Path $projectRoot ("data\campaigns\" + $RunId)
$localPhysicsFile = Join-Path $projectRoot "p53cad\engine\physics_validation.py"

if (-not (Test-Path $localCampaignDir)) {
    throw "Local campaign directory not found: $localCampaignDir"
}
if (-not (Test-Path (Join-Path $localCampaignDir "top30.parquet"))) {
    throw "top30.parquet missing in $localCampaignDir"
}
if (-not (Test-Path $localPhysicsFile)) {
    throw "File missing: $localPhysicsFile"
}

$target = "$RemoteUser@$RemoteHost"
$remoteCampaignBase = "$RemoteProjectDir/data/campaigns"
$remoteCampaignDir = "$remoteCampaignBase/$RunId"
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$remoteLogRel = "logs/karbon_offload_validate_${RunId}_${timestamp}.log"
$remoteLogAbs = "$RemoteProjectDir/$remoteLogRel"
$skipEsmfoldArg = if ($IncludeEsmfold) { "" } else { "--skip-esmfold" }

$sshCommon = @()
if ($KexAlgorithms) {
    $sshCommon += @("-o", "KexAlgorithms=$KexAlgorithms")
}
if ($ConnectTimeoutSec -gt 0) {
    $sshCommon += @("-o", "ConnectTimeout=$ConnectTimeoutSec")
}
if ($BatchMode) {
    $sshCommon += @("-o", "BatchMode=yes")
}
if ($AcceptNewHostKey) {
    $sshCommon += @("-o", "StrictHostKeyChecking=accept-new")
}
if ($IdentityFile) {
    $sshCommon += @("-i", $IdentityFile)
}
if ($ProxyJump) {
    $sshCommon += @("-J", $ProxyJump)
}

Write-Host "Project root        : $projectRoot"
Write-Host "Run ID              : $RunId"
Write-Host "Remote target       : $target"
Write-Host "Remote project dir  : $RemoteProjectDir"
Write-Host "Remote log          : $remoteLogAbs"
Write-Host ""

if (-not $NoStopLocal) {
    Write-Host "[1/6] Stopping local validate processes for run $RunId..."
    Get-CimInstance Win32_Process |
        Where-Object {
            $_.Name -eq "python.exe" -and
            $_.CommandLine -match "p53cad\.cli\.main validate" -and
            $_.CommandLine -match [Regex]::Escape($RunId)
        } |
        ForEach-Object {
            Stop-Process -Id $_.ProcessId -Force
            Write-Host ("  stopped pid " + $_.ProcessId)
        }
} else {
    Write-Host "[1/6] Skipping local process stop (-NoStopLocal set)."
}

Write-Host "[2/6] Ensuring remote directories exist..."
$argsSshMkdir = $sshCommon + @(
    $target,
    "mkdir -p '$RemoteProjectDir' '$remoteCampaignBase' '$RemoteProjectDir/logs' '$RemoteProjectDir/p53cad/engine'"
)
Invoke-Checked -Exe "ssh" -Args $argsSshMkdir

Write-Host "[3/6] Syncing campaign artifacts to remote..."
Invoke-Checked -Exe "scp" -Args ($sshCommon + @(
    "-r",
    $localCampaignDir,
    "${target}:$remoteCampaignBase/"
))

Write-Host "[4/6] Syncing patched physics_validation.py to remote..."
Invoke-Checked -Exe "scp" -Args ($sshCommon + @(
    $localPhysicsFile,
    "${target}:$RemoteProjectDir/p53cad/engine/physics_validation.py"
))

Write-Host "[5/6] Running remote validation on Karbon (CPU/OpenCL path)..."
$remoteCmd = @(
    "cd '$RemoteProjectDir'",
    "export PYTHONIOENCODING=utf-8",
    "export KMP_DUPLICATE_LIB_OK=TRUE",
    "$RemotePython -X utf8 -u -m p53cad.cli.main validate --run-id $RunId --device cpu $skipEsmfoldArg 2>&1 | tee '$remoteLogRel'"
) -join " && "
Invoke-Checked -Exe "ssh" -Args @($sshCommon + @($target, $remoteCmd))

Write-Host "[6/6] Pulling results back to local..."
Invoke-Checked -Exe "scp" -Args ($sshCommon + @(
    "${target}:$remoteCampaignDir/physics_validation.json",
    "$localCampaignDir\"
))

# Optional caches/logs: best effort copy if present
foreach ($entry in @("binding_sim_cache", "esmfold_cache")) {
    try {
        Invoke-Checked -Exe "scp" -Args ($sshCommon + @(
            "-r",
            "${target}:$remoteCampaignDir/$entry",
            "$localCampaignDir\"
        ))
    } catch {
        Write-Warning "Optional pull skipped ($entry): $($_.Exception.Message)"
    }
}

try {
    Invoke-Checked -Exe "scp" -Args ($sshCommon + @(
        "${target}:$remoteLogAbs",
        (Join-Path $projectRoot "logs\")
    ))
} catch {
    Write-Warning "Could not pull remote run log: $($_.Exception.Message)"
}

Write-Host ""
Write-Host "Done. Local result:"
Write-Host ("  " + (Join-Path $localCampaignDir "physics_validation.json"))
Write-Host "Latest local workflow log:"
Write-Host ("  " + (Join-Path $projectRoot "logs\p53cad_workflow.log"))
