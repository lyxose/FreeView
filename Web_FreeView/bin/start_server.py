from __future__ import annotations

import argparse
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import http.client


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

    def do_POST(self):
        # 代理转发 /api/events 请求
        if self.path.startswith("/api/events"):
            content_length = int(self.headers.get('Content-Length', 0))
            post_data = self.rfile.read(content_length)

            try:
                # 使用 http.client 转发请求
                conn = http.client.HTTPConnection("192.168.71.50", 80, timeout=5)
                headers = {"Content-Type": "application/json"}
                conn.request("POST", "/api/events", body=post_data, headers=headers)
                resp = conn.getresponse()

                # 返回响应给客户端
                self.send_response(resp.status)
                for k, v in resp.getheaders():
                    if k.lower() not in ["content-encoding", "transfer-encoding", "connection"]:
                        self.send_header(k, v)
                self.end_headers()
                self.wfile.write(resp.read())
                conn.close()
            except Exception as e:
                body = f"Proxy error: {e}".encode("utf-8")
                self.send_response(502)
                self.send_header("Content-Type", "text/plain; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
        else:
            self.send_response(404)
            self.end_headers()

    def log_message(self, format, *args):  # 保留原有日志功能
        print(f"[INFO] {self.address_string()} - {format % args}")


def main() -> int:
    parser = argparse.ArgumentParser(description='Static file server with Tobii API proxy')
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
