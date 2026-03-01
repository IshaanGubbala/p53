param(
    [string]$RemoteUser = "ishaan",
    [string]$RemoteHost = "karbon800",
    [string]$RemoteProjectDir = "/home/ishaan/p53",
    [string]$RemotePython = "python",
    [switch]$IncludeEsmfold,
    [switch]$AcceptNewHostKey,
    [switch]$NoStopLocal
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Log {
    param([string]$Message)
    $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    "$ts $Message" | Tee-Object -FilePath $script:LogPath -Append
}

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$campaignsDir = Join-Path $projectRoot "data\campaigns"
$logDir = Join-Path $projectRoot "logs"
$statePath = Join-Path $logDir "auto_offload_state.json"
$script:LogPath = Join-Path $logDir "auto_offload.log"
$offloadScript = Join-Path $PSScriptRoot "offload_validate_to_karbon.ps1"

New-Item -ItemType Directory -Path $logDir -Force | Out-Null

if (-not (Test-Path $offloadScript)) {
    throw "Offload script not found: $offloadScript"
}
if (-not (Test-Path $campaignsDir)) {
    Write-Log "No campaigns directory at $campaignsDir; nothing to do."
    exit 0
}

$completedRuns = Get-ChildItem $campaignsDir -Directory |
    Where-Object { Test-Path (Join-Path $_.FullName "top30.parquet") } |
    Sort-Object LastWriteTime -Descending

if ($completedRuns.Count -eq 0) {
    Write-Log "No completed campaign runs found."
    exit 0
}

$latest = $completedRuns[0]
$runId = $latest.Name

$state = @{ processed_run_ids = @() }
if (Test-Path $statePath) {
    try {
        $loaded = Get-Content $statePath -Raw | ConvertFrom-Json
        if ($loaded.processed_run_ids) {
            $state.processed_run_ids = @($loaded.processed_run_ids)
        }
    } catch {
        Write-Log "State file parse failed; resetting state."
    }
}

if ($state.processed_run_ids -contains $runId) {
    Write-Log "Run $runId already processed; skipping."
    exit 0
}

Write-Log "Starting offload for latest run: $runId"

$args = @(
    "-ExecutionPolicy", "Bypass",
    "-File", $offloadScript,
    "-RemoteUser", $RemoteUser,
    "-RemoteHost", $RemoteHost,
    "-RemoteProjectDir", $RemoteProjectDir,
    "-RemotePython", $RemotePython,
    "-RunId", $runId
)
if ($IncludeEsmfold) { $args += "-IncludeEsmfold" }
if ($AcceptNewHostKey) { $args += "-AcceptNewHostKey" }
if ($NoStopLocal) { $args += "-NoStopLocal" }

try {
    & powershell @args
    if ($LASTEXITCODE -ne 0) {
        throw "offload script exited with code $LASTEXITCODE"
    }
    $state.processed_run_ids = @($state.processed_run_ids + $runId | Select-Object -Unique)
    $state | ConvertTo-Json -Depth 4 | Set-Content -Path $statePath -Encoding utf8
    Write-Log "Offload completed for $runId"
    exit 0
} catch {
    Write-Log ("Offload failed for " + $runId + ": " + $_.Exception.Message)
    exit 1
}
