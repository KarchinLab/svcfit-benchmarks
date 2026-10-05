#!/usr/bin/env bash
# Compatibility entrypoint for the project-owned ROCKFISH tool environments.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"

exec "$repo_root/tools/setup_rockfish_tool_envs.sh" "$@"
