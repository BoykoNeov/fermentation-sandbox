"""Start the console the way someone who has never used a terminal would want it started.

``uv run streamlit run app/main.py`` is the honest command and it stays documented, but it
assumes three things a first-time user does not have: that the interface dependencies are
already installed, that port 8501 is free, and that a bare framework error message is a
useful thing to be handed. This script is what the double-clickable wrappers at the repo
root call, and it deals with all three.

Three deliberate choices:

*It installs nothing behind your back.* If ``uv`` is missing the wrappers say so and give the
address to get it; nothing here downloads an installer. If the interface dependencies are
missing this prints the one command that adds them rather than running it — ``uv run --group
ui`` in the wrappers has normally already done it.

*It picks a port that is actually free.* Streamlit refuses to start on a taken port and says
so in one line that scrolls past; on a machine where anything else is already serving on
8501 — including a second copy of this console — that reads as "it did not work". The port
is chosen by asking each candidate whether anything is already answering on it.

*It runs headless and opens the browser itself.* Left to open the browser on its own,
Streamlit first prompts on the terminal for an email address, and a prompt nobody expects
looks like a hang. Headless skips the prompt; the browser is opened here instead.
"""

from __future__ import annotations

import importlib.util
import socket
import subprocess
import sys
import threading
import webbrowser
from pathlib import Path

#: Ports tried in order. 8501 is Streamlit's own default, so a console started this way and
#: one started with the documented command land in the same place when nothing is in the way.
PORT_CANDIDATES: tuple[int, ...] = (8501, 8502, 8503, 8504, 8505, 8506, 8507, 8508)

#: How long to wait for an answer when checking whether a port is already serving.
CONNECT_TIMEOUT_SECONDS = 0.35

#: How long to give the server before pointing a browser at it. The page itself waits for the
#: server, so this only decides whether the first paint is the app or a moment of nothing.
BROWSER_DELAY_SECONDS = 2.5

APP_PATH = Path(__file__).resolve().parent / "main.py"


def say(message: str = "") -> None:
    """Print, and actually put it on the screen.

    Two hazards, both of which only show up away from a developer's terminal.

    Python block-buffers stdout when it is not a terminal, and the server this launches
    writes straight to the same handle. Unflushed, every sentence here arrives *after* the
    framework's own output when the window is logged or piped -- so the explanation of what
    is happening turns up underneath the thing it was explaining.

    And everything said here is ASCII on purpose. A Windows console renders anything, but the
    moment the same output is redirected to a file, Python falls back to the machine's legacy
    encoding and a single em dash raises ``UnicodeEncodeError`` -- a launcher that crashes on
    its own greeting. Measured in this repo, on an arrow in a status line.
    """
    print(message, flush=True)


def port_is_free(port: int, host: str = "127.0.0.1") -> bool:
    """True if nothing is already answering on this port.

    This asks by *connecting*, not by binding, and the difference is not academic. The
    reflex probe — bind the port and see if it raises — was measured against a console
    already serving on 8613 and reported the port free, twice, with and without
    ``SO_REUSEADDR`` on the probe. Windows hands the port over when the *existing* listener
    set that option, which Uvicorn does, so a bind proves nothing about whether anyone is
    there. A refused connection does.
    """
    try:
        with socket.create_connection((host, port), timeout=CONNECT_TIMEOUT_SECONDS):
            return False
    except OSError:
        return True


def choose_port(candidates: tuple[int, ...] = PORT_CANDIDATES) -> int | None:
    """The first free port among the candidates, or ``None`` if every one is taken."""
    return next((port for port in candidates if port_is_free(port)), None)


def missing_interface_dependencies() -> list[str]:
    """Which of the interface packages are not importable from the running interpreter."""
    needed = ("streamlit", "plotly", "pandas")
    return [name for name in needed if importlib.util.find_spec(name) is None]


def command(port: int) -> list[str]:
    """The Streamlit command line, with the choices this launcher makes baked in."""
    return [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(APP_PATH),
        "--server.port",
        str(port),
        # Headless so the framework never stops to ask the terminal for an email address; the
        # browser is opened from here instead.
        "--server.headless",
        "true",
        "--browser.gatherUsageStats",
        "false",
    ]


def main() -> int:
    missing = missing_interface_dependencies()
    if missing:
        say("The console needs a few extra packages that are not installed yet:")
        say("    " + ", ".join(missing))
        say()
        say("Install them with:    uv sync --group ui")
        return 1

    port = choose_port()
    if port is None:
        tried = ", ".join(str(p) for p in PORT_CANDIDATES)
        say("Every port the console normally uses is busy on this machine.")
        say(f"Tried: {tried}.")
        say()
        say("Something else is serving on all of them - very often another copy of this")
        say("console that is still running. Close it, or start this one on a port of your")
        say("own choosing with:")
        say("    uv run streamlit run app/main.py --server.port 8600")
        return 1

    url = f"http://localhost:{port}"
    say()
    say("  The Fermentation Console is starting.")
    say(f"  It will open in your browser at  {url}")
    say()
    say("  Leave this window open while you use it. Closing it, or pressing Ctrl+C here,")
    say("  shuts the console down.")
    say()

    opener = threading.Timer(BROWSER_DELAY_SECONDS, webbrowser.open, args=(url,))
    opener.daemon = True
    opener.start()
    try:
        return subprocess.call(command(port))
    except KeyboardInterrupt:
        say("\n  Console stopped.")
        return 0
    finally:
        opener.cancel()


if __name__ == "__main__":
    raise SystemExit(main())
