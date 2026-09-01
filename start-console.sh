#!/usr/bin/env bash
# Run this file to start the Fermentation Console:  ./start-console.sh
#
# It does three things: check that uv is installed, let uv install what the console needs,
# and hand over to app/start_console.py, which picks a free port and opens the browser.
# Nothing here installs anything system-wide, and nothing is uploaded anywhere.

set -euo pipefail
cd "$(dirname "$0")"

if ! command -v uv >/dev/null 2>&1; then
    cat <<'MSG'

  The console is started by a tool called uv, which is not installed on this machine.

  Install it once from:  https://docs.astral.sh/uv/getting-started/installation/
  Then run this file again.

MSG
    exit 1
fi

cat <<'MSG'

  Getting the console ready. The first time takes a minute or two while the
  pieces are downloaded; after that it is a few seconds.

MSG

exec uv run --group ui python app/start_console.py
