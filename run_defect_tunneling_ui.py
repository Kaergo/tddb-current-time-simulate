#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import functools
import http.server
from pathlib import Path
import socketserver
import webbrowser


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


class ThreadedTCPServer(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch the local defect tunneling visualization UI.")
    parser.add_argument("--host", default="127.0.0.1", help="Bind host.")
    parser.add_argument("--port", type=int, default=8765, help="Preferred port.")
    parser.add_argument("--no-open", action="store_true", help="Do not open the browser automatically.")
    return parser.parse_args()


def serve(host: str, preferred_port: int, no_open: bool) -> None:
    root = Path(__file__).resolve().parent
    html_name = "defect_tunneling_ui.html"
    html_path = root / html_name
    if not html_path.exists():
        raise FileNotFoundError(f"Missing UI file: {html_path}")

    handler = functools.partial(QuietHandler, directory=str(root))
    last_error: OSError | None = None
    for port in range(preferred_port, preferred_port + 20):
        try:
            with ThreadedTCPServer((host, port), handler) as httpd:
                url = f"http://{host}:{port}/{html_name}"
                print(f"Defect tunneling UI: {url}", flush=True)
                print("Press Ctrl+C to stop.", flush=True)
                if not no_open:
                    webbrowser.open(url)
                httpd.serve_forever()
                return
        except OSError as exc:
            last_error = exc
            continue
    raise OSError(f"Could not bind any port from {preferred_port} to {preferred_port + 19}") from last_error


def main() -> None:
    args = parse_args()
    serve(args.host, args.port, args.no_open)


if __name__ == "__main__":
    main()
