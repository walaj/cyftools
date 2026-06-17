#!/usr/bin/env python3
"""Open the cyfview GUI on a CYF table — no clicking, no manual conversion.

    python3 viewer/view.py                 # just serve (cwd as root); load files by URL
    python3 viewer/view.py cells.byf       # serve AND auto-open that file
    python3 viewer/view.py cells.byf 9000   # ... on a specific port

The GUI needs the packed columnar `.cyfv` render format. This server builds one
*under the hood* and caches it in a temp dir keyed on the source file, so you
never see or manage a `.cyfv`. Two ways to load:

  • auto-open:  pass a `.byf` (or `.cyf`/`.ocyf`/`.cyfv`) on the command line.
  • by URL:     request any table path and it is converted+cached on the fly, e.g.
        http://localhost:8765/viewer/cyfview.html?file=/test/LSP10353.byf
    Paths in `?file=` resolve under the server root (the cwd you launched from,
    or $CYFVIEW_ROOT). Just swap the path to load a different file — no restart.

Edit the `.byf` (e.g. with `cyftools pheno`) and refresh: the cache key changes,
so the pack rebuilds on demand. A `.cyfv` is served as-is (no conversion).

Server mode (containers / Cloud Run): set $PORT and it binds 0.0.0.0:$PORT and
serves without opening a browser — users upload a `.byf` in the page itself and
it is converted server-side (POST /convert), so they need no local cyftools.

Env: PORT (server mode), CYFTOOLS=/path/to/cyftools, CYFVIEW_CACHE=/cachedir,
CYFVIEW_ROOT=/served/root. Ctrl-C stops.
"""
import os, sys, time, hashlib, tempfile, shutil, subprocess
import json, uuid, datetime
import http.server, socketserver, webbrowser, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HTML = os.path.join(HERE, "cyfview.html")

TABLE_EXTS = (".byf", ".cyf", ".ocyf")     # source tables -> need conversion
PACK_EXTS  = TABLE_EXTS + (".cyfv",)        # anything the loader may request
CACHE_FMT  = "v5"                            # bump if the .cyfv layout changes
CACHE_DIR  = os.environ.get("CYFVIEW_CACHE") or os.path.join(tempfile.gettempdir(), "cyfview-cache")
CACHE_TTL  = 7 * 24 * 3600                   # prune cached packs older than this
ROOT       = os.path.abspath(os.environ.get("CYFVIEW_ROOT") or os.getcwd())

# Big-file path. Uploads go to an object store (a GCS bucket in the cloud, or the
# local cache dir for dev) instead of through the server, so multi-GB files are not
# bound by any request-body limit. MAX_CELLS caps what the browser has to render:
# beyond it the pack is built from a streaming `subsample` (full-res stays in the
# store), keeping the downloaded pack and the in-browser buffers bounded.
BUCKET     = os.environ.get("CYFVIEW_BUCKET")          # set -> GCS store; unset -> local
MAX_CELLS  = int(os.environ.get("CYFVIEW_MAX_CELLS", "4000000"))
SIGN_TTL   = datetime.timedelta(minutes=30)            # signed-URL lifetime
STORE      = None                                       # set in main()


def find_cyftools():
    """Locate the cyftools binary: $CYFTOOLS, then build/, then PATH."""
    env = os.environ.get("CYFTOOLS")
    if env and os.path.isfile(env) and os.access(env, os.X_OK):
        return env
    for c in (os.path.join(REPO, "build", "cyftools"), os.path.join(REPO, "cyftools")):
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return shutil.which("cyftools")


def cache_path_for(src, exporter):
    """Deterministic cache path for this source's current state.

    The key folds in the source path, its mtime+size, the exporter's mtime, and
    the format tag — so editing the .byf or rebuilding cyftools yields a new key
    (and a fresh pack), while an unchanged file reuses the existing one.
    """
    st = os.stat(src)
    ex_mtime = os.stat(exporter).st_mtime_ns if exporter and os.path.isfile(exporter) else 0
    key = "\0".join([os.path.abspath(src), str(st.st_mtime_ns), str(st.st_size),
                     str(ex_mtime), CACHE_FMT])
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    base = os.path.splitext(os.path.basename(src))[0]
    return os.path.join(CACHE_DIR, f"{base}.{h}.cyfv")


def prune_cache():
    try:
        now = time.time()
        for f in os.listdir(CACHE_DIR):
            p = os.path.join(CACHE_DIR, f)
            if f.endswith(".cyfv") and now - os.path.getmtime(p) > CACHE_TTL:
                os.remove(p)
    except OSError:
        pass


def count_cells(src, exporter):
    """Total cells in `src` (streaming, cheap). -1 if it can't be determined."""
    try:
        r = subprocess.run([exporter, "count", src], capture_output=True, text=True, check=True)
        return int(r.stdout.strip().split()[0])
    except (subprocess.CalledProcessError, ValueError, IndexError):
        return -1


def build_pack(src, dest, exporter):
    """Export `src` to a CYFV pack at `dest`, capping the rendered cells at MAX_CELLS.

    Over the cap, the pack is built from a streaming `subsample -r` (RAM-flat, no full
    load) so the browser payload stays bounded no matter how large the input is.
    Returns (total_cells, rendered_cells, subsampled). Publishes `dest` atomically.
    """
    if not exporter:
        raise RuntimeError("cyftools not found — set $CYFTOOLS or build it (build/cyftools)")
    total = count_cells(src, exporter)
    tmp = f"{dest}.{os.getpid()}.tmp"
    try:
        if total > MAX_CELLS:
            rate = MAX_CELLS / total
            sub = subprocess.Popen([exporter, "subsample", src, "-", "-r", f"{rate:.6f}"],
                                   stdout=subprocess.PIPE)
            exp = subprocess.Popen([exporter, "export", "-", tmp], stdin=sub.stdout)
            sub.stdout.close()                       # let export own the read end
            exp.communicate()
            sub.wait()
            if sub.returncode or exp.returncode:
                raise subprocess.CalledProcessError(exp.returncode or sub.returncode, "subsample|export")
            rendered, subsampled = int(total * rate), True
        else:
            subprocess.run([exporter, "export", src, tmp], check=True)
            rendered, subsampled = total, False
    except subprocess.CalledProcessError as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise RuntimeError(f"cyftools pack build failed ({e.returncode})")
    os.replace(tmp, dest)                            # atomic publish
    return total, rendered, subsampled


def ensure_pack(src, exporter):
    """Return a path to a current .cyfv for `src`, building into cache if needed."""
    if src.lower().endswith(".cyfv"):
        return src                                   # already a pack: serve as-is
    if not exporter:
        raise RuntimeError("cyftools not found — set $CYFTOOLS or build it (build/cyftools)")
    os.makedirs(CACHE_DIR, exist_ok=True)
    out = cache_path_for(src, exporter)
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        return out                                   # cache hit
    sys.stderr.write(f"cyfview: converting {os.path.basename(src)} -> cache ...\n")
    build_pack(src, out, exporter)                   # applies the render cap
    prune_cache()
    return out


def resolve_source(reqpath):
    """Map a requested URL path to a real table/pack file under ROOT (or absolute)."""
    rel = urllib.parse.unquote(reqpath).lstrip("/")
    for cand in (os.path.normpath(os.path.join(ROOT, rel)),
                 os.path.normpath(urllib.parse.unquote(reqpath))):
        if os.path.isfile(cand) and cand.lower().endswith(PACK_EXTS):
            return cand
    return None


class LocalStore:
    """Dev object store: blobs live under CACHE_DIR; the browser PUTs uploads to
    /blob/uploads/<id> on this same server and fetches packs from /blob/packs/<id>."""
    kind = "local"

    def upload_path(self, blob_id):
        return os.path.join(CACHE_DIR, "uploads", blob_id)

    def upload_target(self, blob_id):
        return {"url": f"/blob/uploads/{blob_id}", "method": "PUT", "headers": {}}

    def fetch_to_local(self, blob_id):
        return self.upload_path(blob_id)             # already local

    def publish_pack(self, blob_id, local_pack):
        dst = os.path.join(CACHE_DIR, "packs", f"{blob_id}.cyfv")
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.abspath(local_pack) != os.path.abspath(dst):
            os.replace(local_pack, dst)

    def pack_url(self, blob_id):
        return f"/blob/packs/{blob_id}.cyfv"


class GcsStore:
    """Cloud object store: the browser uploads straight to gs://BUCKET/uploads/<id>
    via a signed PUT URL (no request-size limit, never through the server); packs land
    in gs://BUCKET/packs/<id>.cyfv and are fetched via a signed GET URL. The server only
    streams blobs to/from local disk so cyftools can read/write real files."""
    kind = "gcs"

    def __init__(self, bucket):
        from google.cloud import storage            # lazy: only needed in cloud mode
        import google.auth
        from google.auth.transport import requests as ga_requests
        self._bucket = storage.Client().bucket(bucket)
        # Credentials for V4 signing. A mounted service-account key signs locally; on
        # Cloud Run's keyless default SA there is no private key, so we sign via the IAM
        # API (needs roles/iam.serviceAccountTokenCreator on the SA itself).
        self._creds, _ = google.auth.default()
        self._auth_req = ga_requests.Request()
        self._sa_email = self._lookup_sa_email()

    def _lookup_sa_email(self):
        sa = getattr(self._creds, "service_account_email", None)
        if sa and sa != "default":
            return sa
        try:                                         # GCE/Cloud Run metadata server
            req = urllib.request.Request(
                "http://metadata.google.internal/computeMetadata/v1/"
                "instance/service-accounts/default/email",
                headers={"Metadata-Flavor": "Google"})
            return urllib.request.urlopen(req, timeout=2).read().decode().strip()
        except Exception:
            return None

    def _signed(self, blob, method, content_type=None):
        kw = dict(version="v4", expiration=SIGN_TTL, method=method)
        if content_type:
            kw["content_type"] = content_type
        try:
            return blob.generate_signed_url(**kw)    # local key path
        except Exception:                            # keyless: sign through IAM signBlob
            self._creds.refresh(self._auth_req)
            return blob.generate_signed_url(service_account_email=self._sa_email,
                                            access_token=self._creds.token, **kw)

    def upload_target(self, blob_id):
        url = self._signed(self._bucket.blob(f"uploads/{blob_id}"),
                           "PUT", content_type="application/octet-stream")
        return {"url": url, "method": "PUT", "headers": {"Content-Type": "application/octet-stream"}}

    def fetch_to_local(self, blob_id):
        dst = os.path.join(CACHE_DIR, "uploads", blob_id)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        self._bucket.blob(f"uploads/{blob_id}").download_to_filename(dst)
        return dst

    def publish_pack(self, blob_id, local_pack):
        self._bucket.blob(f"packs/{blob_id}.cyfv").upload_from_filename(
            local_pack, content_type="application/octet-stream")
        for p in (local_pack, os.path.join(CACHE_DIR, "uploads", blob_id)):
            try: os.remove(p)                        # ephemeral: reclaim temp disk
            except OSError: pass

    def pack_url(self, blob_id):
        return self._signed(self._bucket.blob(f"packs/{blob_id}.cyfv"), "GET")


def make_store():
    if BUCKET:
        try:
            return GcsStore(BUCKET)
        except ImportError:
            sys.exit("CYFVIEW_BUCKET is set but google-cloud-storage is not installed "
                     "(pip install google-cloud-storage)")
    return LocalStore()


def main():
    global STORE
    args = sys.argv[1:]
    launch_src, port = None, 8765
    for a in args:
        if a.isdigit():
            port = int(a)
        else:
            launch_src = os.path.abspath(a)
    # Server mode: when $PORT is set (Cloud Run / containers) bind 0.0.0.0:$PORT and
    # never open a browser. Otherwise it is a local launch on 127.0.0.1.
    server_mode = bool(os.environ.get("PORT"))
    if server_mode:
        port = int(os.environ["PORT"])
    host = "0.0.0.0" if server_mode else "127.0.0.1"
    STORE = make_store()
    if not os.path.isfile(HTML):
        sys.exit(f"cyfview.html not found next to this script ({HTML})")
    exporter = find_cyftools()
    if launch_src:
        if not os.path.isfile(launch_src):
            sys.exit(f"no such file: {launch_src}")
        if not launch_src.lower().endswith(PACK_EXTS):
            sys.exit(f"unsupported input (want .byf/.cyf/.ocyf/.cyfv): {launch_src}")
        try:
            ensure_pack(launch_src, exporter)        # build cache up front
        except RuntimeError as e:
            sys.exit(str(e))

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass
        def _send_file(self, path, ctype):
            try:
                with open(path, "rb") as f:
                    body = f.read()
            except OSError:
                self.send_error(404); return
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(body)
        def _serve_pack(self, src):
            try:
                pack = ensure_pack(src, exporter)
            except RuntimeError as e:
                self.send_error(500, str(e)); return
            self._send_file(pack, "application/octet-stream")
        def _json_body(self):
            n = int(self.headers.get("Content-Length", 0) or 0)
            raw = self.rfile.read(n) if n else b""
            return json.loads(raw or b"{}")
        def _send_json(self, obj, code=200):
            body = json.dumps(obj).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(body)
        def _serve_blob(self, p):                    # GET /blob/<rel> from CACHE_DIR (local store)
            rel = urllib.parse.unquote(p[len("/blob/"):])
            full = os.path.normpath(os.path.join(CACHE_DIR, rel))
            if not full.startswith(os.path.abspath(CACHE_DIR) + os.sep):
                self.send_error(403); return
            self._send_file(full, "application/octet-stream")
        def do_GET(self):
            p = urllib.parse.urlparse(self.path).path
            if p in ("/", "/cyfview.html", "/viewer/cyfview.html"):
                self._send_file(HTML, "text/html")
            elif p == "/data.byf" and launch_src:         # the file passed on the CLI
                self._serve_pack(launch_src)
            elif STORE.kind == "local" and p.startswith("/blob/"):   # local object store
                self._serve_blob(p)
            elif p.lower().endswith(PACK_EXTS):           # any table/pack requested by URL
                src = resolve_source(p)
                if src:
                    self._serve_pack(src)
                else:
                    self.send_error(404, "no such table under server root")
            else:
                self.send_error(404)

        def do_PUT(self):
            # local object store: the browser PUTs a .byf here (streamed to disk).
            p = urllib.parse.urlparse(self.path).path
            if STORE.kind != "local" or not p.startswith("/blob/uploads/"):
                self.send_error(404); return
            blob_id = p[len("/blob/uploads/"):]
            if not blob_id or not blob_id.isalnum():
                self.send_error(400, "bad blob id"); return
            dst = STORE.upload_path(blob_id)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            remaining = int(self.headers.get("Content-Length", 0) or 0)
            try:
                with open(dst, "wb") as f:
                    while remaining > 0:
                        chunk = self.rfile.read(min(1 << 20, remaining))
                        if not chunk:
                            break
                        f.write(chunk); remaining -= len(chunk)
            except Exception:
                self.send_error(500, "write failed"); return
            self.send_response(200)
            self.send_header("Content-Length", "0")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()

        def _api_upload_url(self):
            # hand the browser a place to PUT a .byf (signed GCS URL, or local /blob)
            try:
                self._json_body()                    # parse/validate (name,size ignored)
            except Exception:
                self.send_error(400, "bad json"); return
            blob_id = uuid.uuid4().hex
            tgt = STORE.upload_target(blob_id)
            self._send_json({"id": blob_id, "put_url": tgt["url"],
                             "method": tgt["method"], "headers": tgt["headers"]})

        def _api_pack(self):
            # convert an already-uploaded blob to a (render-capped) pack in the store
            if not exporter:
                self.send_error(500, "cyftools not found"); return
            try:
                blob_id = self._json_body().get("id", "")
            except Exception:
                self.send_error(400, "bad json"); return
            if not blob_id or not blob_id.isalnum():
                self.send_error(400, "missing/bad id"); return
            try:
                src = STORE.fetch_to_local(blob_id)
                if not (os.path.isfile(src) and os.path.getsize(src) > 0):
                    self.send_error(404, "upload not found"); return
                os.makedirs(CACHE_DIR, exist_ok=True)
                tmp = os.path.join(CACHE_DIR, f"pack.{blob_id}.cyfv")
                total, rendered, subsampled = build_pack(src, tmp, exporter)
                STORE.publish_pack(blob_id, tmp)
            except RuntimeError as e:
                self.send_error(500, str(e)); return
            self._send_json({"pack_url": STORE.pack_url(blob_id), "cells": total,
                             "rendered": rendered, "subsampled": subsampled})

        def do_POST(self):
            p = urllib.parse.urlparse(self.path).path
            if p == "/api/upload-url":
                return self._api_upload_url()
            if p == "/api/pack":
                return self._api_pack()
            # fallback: small-file direct convert (browser POSTs raw .byf bytes)
            if p != "/convert":
                self.send_error(404); return
            if not exporter:
                self.send_error(500, "cyftools not found"); return
            try:
                body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
            except Exception:
                self.send_error(400, "bad upload"); return
            if not body:
                self.send_error(400, "empty upload"); return
            os.makedirs(CACHE_DIR, exist_ok=True)
            h = hashlib.sha1(body).hexdigest()[:16]
            out = os.path.join(CACHE_DIR, f"upload.{h}.cyfv")
            if not (os.path.isfile(out) and os.path.getsize(out) > 0):
                src = os.path.join(CACHE_DIR, f"upload.{h}.byf")
                with open(src, "wb") as fh:
                    fh.write(body)
                try:
                    subprocess.run([exporter, "export", src, out], check=True)
                except subprocess.CalledProcessError as e:
                    try: os.remove(out)
                    except OSError: pass
                    self.send_error(500, f"cyftools export failed ({e.returncode})"); return
                finally:
                    try: os.remove(src)
                    except OSError: pass
            self._send_file(out, "application/octet-stream")

    http.server.ThreadingHTTPServer.allow_reuse_address = True
    if server_mode:
        httpd = http.server.ThreadingHTTPServer((host, port), Handler)   # fixed $PORT
        tryport = port
        print(f"cyfview: serving on {host}:{port} (root: {ROOT})", flush=True)
        print("  upload a .byf in the browser; it is converted server-side", flush=True)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped")
        return

    for tryport in range(port, port + 20):
        try:
            httpd = http.server.ThreadingHTTPServer((host, tryport), Handler)
            break
        except OSError:
            continue
    else:
        sys.exit("no free port found")

    base = f"http://localhost:{tryport}"
    print(f"cyfview: serving (root: {ROOT})")
    if launch_src:
        url = f"{base}/?file=/data.byf"
        print(f"  opening {os.path.basename(launch_src)}: {url}")
    else:
        url = f"{base}/viewer/cyfview.html?file=/path/to/cells.byf"
        print(f"  load by URL: {url}")
    print("  (edit the .byf and refresh to rebuild · Ctrl-C to stop)")
    if launch_src:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")


if __name__ == "__main__":
    main()
