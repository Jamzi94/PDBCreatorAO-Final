#!/usr/bin/env bash
# usage: source ./setup_env.sh
set -euo pipefail
ENV_NAME="myproj"
YAML_FILE="environment.yml"

command_exists() {
  command -v "$1" &>/dev/null
}

# Locate or install micromamba if conda not found
if command_exists conda; then
  PM=conda
elif command_exists micromamba; then
  PM=micromamba
else
  echo "🔍 No conda found – installing micromamba..."
  MAMBA_DIR="$HOME/.local/micromamba"
  mkdir -p "$MAMBA_DIR"
  OS=$(uname | tr '[:upper:]' '[:lower:]')
  ARCH=$(uname -m)
  [[ $ARCH == x86_64 ]] && ARCH=64 || ARCH=32
  curl -Ls "https://micro.mamba.pm/api/micromamba/${OS}-${ARCH}/latest" \
       | tar -C "$MAMBA_DIR" -xj bin/micromamba
  export PATH="$MAMBA_DIR/bin:$PATH"
  PM=micromamba
fi
echo "🛠️  Using package manager: $PM"

# Initialize shell
eval "$($PM shell hook --shell bash)"

# Remove existing env if present and create new
$PM env remove -y -n "$ENV_NAME" &>/dev/null || true
$PM env create -f "$YAML_FILE" -n "$ENV_NAME"

# Activate in current shell
$PM activate "$ENV_NAME"
echo "✅ Environment '$ENV_NAME' is active. Python: $(which python)"
