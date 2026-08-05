"""Entry point for the KEGGaNOG desktop web application.

Orchestrates the lifecycle of the Uvicorn ASGI server and coordinates
the automated browser launch sequence.
"""

import threading
import time
import webbrowser
from typing import Final

import uvicorn

# Constants for server binding
HOST: Final[str] = "127.0.0.1"
PORT: Final[int] = 8000
BROWSER_DELAY: Final[float] = 1.5


def _open_browser_delayed(url: str, delay: float) -> None:
    """Task executed in a daemon thread to launch the UI after a delay."""
    time.sleep(delay)
    webbrowser.open(url)


def launch() -> None:
    """Launch the FastAPI server and trigger the browser interface."""
    url = f"http://{HOST}:{PORT}"

    # Start browser-opener in a background daemon thread
    browser_thread = threading.Thread(
        target=_open_browser_delayed,
        args=(url, BROWSER_DELAY),
        daemon=True,
    )
    browser_thread.start()

    print(f"\n  {'=' * 40}")
    print(f"  KEGGaNOG web UI is running at: {url}")
    print("  Press Ctrl+C to terminate the server.")
    print(f"  {'=' * 40}\n")

    # Start the blocking ASGI server
    try:
        uvicorn.run(
            "kegganog.app:app",
            host=HOST,
            port=PORT,
            log_level="warning",
            reload=False,
        )
    except KeyboardInterrupt:
        print("\n  Server stopped by user. Goodbye!")
    except Exception as e:
        print(f"\n  Critical server error: {e}")


if __name__ == "__main__":
    launch()
