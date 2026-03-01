param(
    [string]$TaskName = "p53CAD-AutoOffload",
    [int]$EveryMinutes = 15,
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

if ($EveryMinutes -lt 1) {
    throw "EveryMinutes must be >= 1"
}

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$worker = Join-Path $PSScriptRoot "auto_offload_once.ps1"
if (-not (Test-Path $worker)) {
    throw "Worker script not found: $worker"
}

$argList = @(
    "-NoProfile",
    "-ExecutionPolicy", "Bypass",
    "-File", "`"$worker`"",
    "-RemoteUser", $RemoteUser,
    "-RemoteHost", $RemoteHost,
    "-RemoteProjectDir", $RemoteProjectDir,
    "-RemotePython", $RemotePython
)
if ($IncludeEsmfold) { $argList += "-IncludeEsmfold" }
if ($AcceptNewHostKey) { $argList += "-AcceptNewHostKey" }
if ($NoStopLocal) { $argList += "-NoStopLocal" }

$taskCmd = "powershell.exe " + ($argList -join " ")

# Replace existing task if present (ignore if absent).
schtasks /Delete /TN $TaskName /F | Out-Null

schtasks /Create `
    /TN $TaskName `
    /SC MINUTE `
    /MO $EveryMinutes `
    /TR $taskCmd `
    /RL LIMITED `
    /F | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Failed to create scheduled task '$TaskName' (exit $LASTEXITCODE)."
}

Write-Host "Installed scheduled task: $TaskName"
Write-Host "Runs every $EveryMinutes minute(s)."
Write-Host "Command: $taskCmd"
Write-Host ""
Write-Host "Useful commands:"
Write-Host "  schtasks /Query /TN $TaskName /V /FO LIST"
Write-Host "  schtasks /Run /TN $TaskName"
Write-Host "  schtasks /Delete /TN $TaskName /F"
Write-Host "  Get-Content `"$projectRoot\logs\auto_offload.log`" -Wait"
