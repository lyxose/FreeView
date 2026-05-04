from __future__ import annotations

import argparse
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


class ShutdownAwareHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.split('?', 1)[0] == '/__shutdown__':
            body = b'shutting down'
            self.send_response(200)
            self.send_header('Content-Type', 'text/plain; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            self.wfile.flush()
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        return super().do_GET()

    def log_message(self, format, *args):  # noqa: A003 - keep the standard handler signature
        print(f"[INFO] {self.address_string()} - {format % args}")


def main() -> int:
    parser = argparse.ArgumentParser(description='Static file server for FreeView Web Experiment')
    parser.add_argument('--port', type=int, default=8000)
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    handler = partial(ShutdownAwareHandler, directory=str(base_dir))
    server = ThreadingHTTPServer(('127.0.0.1', args.port), handler)
    server.allow_reuse_address = True

    print(f'[INFO] serving {base_dir} on http://127.0.0.1:{args.port}/')
    try:
        server.serve_forever(poll_interval=0.2)
    except KeyboardInterrupt:
        print('[INFO] interrupted by user')
    finally:
        server.server_close()
        print('[INFO] server closed')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())