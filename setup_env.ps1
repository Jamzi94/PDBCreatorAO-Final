<#
  usage: . .\setup_env.ps1   (dot-source to keep environment active)
#>
param(
  [string]$EnvName = "myproj",
  [string]$YamlFile = "environment.yml"
)

function Command-Exists($cmd) {
  Get-Command $cmd -ErrorAction SilentlyContinue | Out-Null
}

# Determine whether to use conda or micromamba
if (-not (Command-Exists conda)) {
  if (-not (Command-Exists micromamba)) {
    Write-Host "🔍 No conda found – installing micromamba..."
    # PowerShell 5.1 doesn't support the C-style ternary operator; use explicit logic
    $osArch = (Get-CimInstance Win32_OperatingSystem).OSArchitecture
    if ($osArch -match "64") { $arch = "64" } else { $arch = "32" }
    $url  = "https://micro.mamba.pm/api/micromamba/win-$arch/latest"
    $out  = Join-Path $env:TEMP 'micromamba.tar.bz2'
    Invoke-WebRequest $url -OutFile $out
    $dest = Join-Path $HOME 'micromamba'
    New-Item -ItemType Directory -Path $dest -Force | Out-Null

    # Extract archive into a temp folder then locate micromamba.exe
    $extractDir = Join-Path $env:TEMP 'micromamba_extract'
    Remove-Item -Recurse -Force $extractDir -ErrorAction SilentlyContinue
    New-Item -ItemType Directory -Path $extractDir -Force | Out-Null
    try {
      tar -xjf $out -C $extractDir 2>$null
    } catch {
      Write-Warning "tar extraction reported an error: $_"
    }

    $exe = Get-ChildItem -Path $extractDir -Recurse -Filter 'micromamba.exe' -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($null -ne $exe) {
      Copy-Item -Path $exe.FullName -Destination (Join-Path $dest 'micromamba.exe') -Force
    } else {
      Write-Host 'micromamba executable not found in archive; attempting direct download...'
      $directUrl = "https://micro.mamba.pm/api/micromamba/win-$arch/latest/micromamba.exe"
      try {
        Invoke-WebRequest $directUrl -OutFile (Join-Path $dest 'micromamba.exe') -ErrorAction Stop
      } catch {
        Write-Error "Failed to obtain micromamba.exe from archive or direct URL. Please install conda/micromamba manually. $_"
        return
      }
    }

    # Ensure the destination is on PATH for this session
    $env:PATH = (Join-Path $dest '') + ";" + $env:PATH
  }
  function conda { micromamba @Args }
}

# Initialize shell hook for either conda or micromamba
$condaCmd = Get-Command conda -ErrorAction SilentlyContinue
$micromambaCmd = Get-Command micromamba -ErrorAction SilentlyContinue
$mmPath = $null
if ($micromambaCmd) { $mmPath = $micromambaCmd.Source }
elseif (Test-Path (Join-Path $dest 'micromamba.exe')) { $mmPath = (Join-Path $dest 'micromamba.exe') }

if ($condaCmd) {
  # If conda is a function we've created that forwards to micromamba, calling
  # "conda 'shell.powershell' 'hook'" will forward the args in the wrong order
  # (micromamba expects: micromamba shell hook -s powershell). Detect function
  # and call the micromamba executable directly when possible.
  $cmdInfo = Get-Command conda -ErrorAction SilentlyContinue
  if ($cmdInfo -and $cmdInfo.CommandType -eq 'Function') {
    $mmExec = Get-Command micromamba -ErrorAction SilentlyContinue
    if ($mmExec) {
      & $mmExec.Source shell hook -s powershell | Out-String | Invoke-Expression
      $manager = $mmExec.Source
    } elseif (Test-Path (Join-Path $dest 'micromamba.exe')) {
      & (Join-Path $dest 'micromamba.exe') shell hook -s powershell | Out-String | Invoke-Expression
      $manager = (Join-Path $dest 'micromamba.exe')
    } else {
      Write-Warning 'Detected conda function but micromamba executable not found; skipping shell hook.'
      $manager = 'conda'
    }
  } else {
    & conda 'shell.powershell' 'hook' | Out-String | Invoke-Expression
    $manager = 'conda'
  }
} elseif ($mmPath) {
  # micromamba expects: micromamba.exe shell hook -s powershell
  & $mmPath shell hook -s powershell | Out-String | Invoke-Expression
  $manager = $mmPath
  # Also attempt to initialize micromamba persistently for future shells
  try {
    & $mmPath shell init --shell powershell --root-prefix=$env:USERPROFILE\/.local/share/mamba | Out-String | Out-Null
    Write-Host "micromamba shell init executed (persistent initialization)."
  } catch {
    Write-Warning "micromamba persistent init failed: $_"
  }
} else {
  Write-Warning 'No conda or micromamba available to initialize shell. Environment commands may fail.'
  $manager = $null
}

# Recreate environment using the detected manager
if ($manager -eq 'conda') {
  conda env remove -y -n $EnvName 2>$null
  conda env create -f $YamlFile -n $EnvName
  conda activate $EnvName
} elseif ($manager) {
  & $manager env remove -y -n $EnvName 2>$null
  & $manager env create -f $YamlFile -n $EnvName
  # Try to activate using shell-initialized function if available
  try {
    & $manager activate $EnvName
  } catch {
    Write-Host "To use the environment, run: `$ $manager activate $EnvName` or follow micromamba's instructions to initialize your shell."
  }
} else {
  Write-Error 'Unable to create environment: no manager available.'
}

Write-Host "✅ Environment '$EnvName' setup attempted. Python: $(Get-Command python -ErrorAction SilentlyContinue).Source"
