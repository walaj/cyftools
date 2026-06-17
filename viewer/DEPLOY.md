# Deploying the cyftools cloud viewer

One Docker image builds the real `cyftools` C++ binary and serves `cyfview.html` plus
`viewer/view.py`. A browser uploads a `.byf`; the server-side binary converts it; the
page renders the result. **There is no second (JavaScript) codec** — every conversion
runs the same `cyftools` you run at the terminal (verified: a container-built pack is
byte-for-byte identical to a local `cyftools export`).

## Two upload paths (the server supports both)

```
small file:   browser ──POST /convert (raw bytes)──▶ view.py ─▶ cyftools ─▶ pack ─▶ browser

big file:     browser ─POST /api/upload-url▶ view.py            (1) get a place to PUT
              browser ─────PUT file─────────▶ object store      (2) upload straight to the bucket
              browser ─POST /api/pack───────▶ view.py           (3) ask for a render pack
                                   view.py ─▶ store→cyftools→store
              browser ◀──signed GET url────── object store      (4) fetch the bounded pack
```

The big-file path is what makes **multi-GB** files work: the browser PUTs the file
**straight into a GCS bucket**, never through the server, so Cloud Run's request-body
limit (~32 MiB) is irrelevant. The server only streams blobs between the bucket and local
disk for `cyftools`. The browser tries the staged path first and falls back to `/convert`
automatically if the server doesn't expose it.

### The render cap (why a few-GB file still renders)

A few-GB `.byf` is tens of millions of cells, and the pack is ~2× the input — too big for
a browser to hold in one buffer or draw smoothly. So the server caps what the browser
gets: it runs a streaming `cyftools count`, and if the table exceeds **`CYFVIEW_MAX_CELLS`**
(default 4,000,000) it builds the pack from `cyftools subsample -r <rate>` (RAM-flat, no
full load). Full resolution stays in the bucket; only the overview is bounded. The page
shows e.g. *"showing 100,000 of 1,620,375 cells (subsampled to fit)"*.

## Prerequisite (once)

The image COPYs `external/` in but not `.git`, so check out submodules first:

```sh
git submodule update --init --recursive
```

## Run locally (Docker)

```sh
docker build -t cyftools-viewer .
docker run --rm -p 8080:8080 cyftools-viewer       # open http://localhost:8080/
```

With no `CYFVIEW_BUCKET`, the server uses a **local object store** (the cache dir): the
browser PUTs to `/blob/uploads/<id>` on the same server, and the whole staged flow works
without any cloud — handy for testing the exact path the cloud uses. There's no
request-size limit on your own host, so any `.byf` works.

## Deploy to Cloud Run with a GCS bucket

### 1. Bucket + permissions

```sh
gcloud storage buckets create gs://MY-CYF-BUCKET --location=us-central1

# the Cloud Run service account needs to read/write objects and to sign URLs:
SA="cyftools-viewer@MY-PROJECT.iam.gserviceaccount.com"
gcloud storage buckets add-iam-policy-binding gs://MY-CYF-BUCKET \
  --member="serviceAccount:$SA" --role="roles/storage.objectAdmin"
gcloud iam service-accounts add-iam-policy-binding "$SA" \
  --member="serviceAccount:$SA" --role="roles/iam.serviceAccountTokenCreator"   # V4 signBlob
```

### 2. Bucket CORS (the browser PUTs/GETs the bucket directly)

```sh
cat > /tmp/cors.json <<'EOF'
[{"origin":["https://cyftools-viewer-XXXX.run.app"],
  "method":["GET","PUT"],
  "responseHeader":["Content-Type"],
  "maxAgeSeconds":3600}]
EOF
gcloud storage buckets update gs://MY-CYF-BUCKET --cors-file=/tmp/cors.json
```

(You'll know the `run.app` origin after the first deploy; update CORS then.)

### 3. Deploy

```sh
gcloud run deploy cyftools-viewer \
  --source . \
  --region us-central1 \
  --service-account "$SA" \
  --allow-unauthenticated \
  --set-env-vars CYFVIEW_BUCKET=MY-CYF-BUCKET,CYFVIEW_MAX_CELLS=4000000 \
  --memory 4Gi --cpu 2 --timeout 900
```

### Sizing / disk

- The server downloads each upload to local disk to run `cyftools`. On Cloud Run that
  disk is **in-memory (tmpfs)** by default, so a multi-GB temp file counts against
  `--memory`. For genuinely large files, mount the bucket as a volume so temps live on
  the bucket, not RAM:
  ```sh
  --add-volume name=data,type=cloud-storage,bucket=MY-CYF-BUCKET \
  --add-volume-mount volume=data,mount-path=/mnt/data \
  --set-env-vars CYFVIEW_CACHE=/mnt/data/cache,CYFVIEW_BUCKET=MY-CYF-BUCKET
  ```
- `--timeout 900` (15 min) covers conversion of large tables; raise toward the 3600 max
  if needed. The signed-URL upload itself isn't bound by this (it's browser→bucket).

### Auth

Drop `--allow-unauthenticated` to require IAM / Google sign-in, or front it with
Identity-Aware Proxy.

## Server behavior / env

| Mode | Trigger | Binds | Browser | Store |
|---|---|---|---|---|
| local launch | `python3 viewer/view.py cells.byf` | `127.0.0.1`, auto port | opens | local |
| server / cloud | `$PORT` set | `0.0.0.0:$PORT` | never | local, or GCS if `CYFVIEW_BUCKET` set |

Env: `CYFTOOLS` (binary path), `CYFVIEW_BUCKET` (GCS bucket → cloud store; unset → local),
`CYFVIEW_MAX_CELLS` (render cap, default 4,000,000), `CYFVIEW_CACHE` (scratch/cache dir),
`CYFVIEW_ROOT` (root for `?file=` path loads).

## Next slices (when you want them)

1. **Generic op endpoint** — `POST /api/run` with a pipeline (`pheno`, `flagroi`,
   `cohort`, …) over a stored blob id, so the whole CLI surface is drivable from the page.
2. **Async jobs** — submit → poll → fetch, for conversions that exceed a request timeout.
3. **Tiled / level-of-detail rendering** — stream cells by viewport instead of one capped
   pack, for interactive detail on the very largest slides.
