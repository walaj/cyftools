# cyftools cloud viewer — one image that builds the C++ binary and serves
# cyfview.html plus a server-side  .byf -> pack  converter (viewer/view.py).
# Users hit it with a browser, upload a .byf, and get a rendered view; the real
# cyftools binary does every conversion, so there is no second (JS) codec.
#
#   local:  docker build -t cyftools-viewer .
#           docker run --rm -p 8080:8080 cyftools-viewer      # http://localhost:8080/
#   cloud:  gcloud run deploy cyftools-viewer --source . --allow-unauthenticated \
#                  --region us-central1 --memory 4Gi --cpu 2 --timeout 600
#
# Prereq: check out the submodules first — external/ is COPYed in, .git is not:
#           git submodule update --init --recursive

# ---------- stage 1: build cyftools (minimal core: no TIFF/Cairo/LDA/OpenMP) ----------
FROM ubuntu:24.04 AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake zlib1g-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
# only the build inputs, so viewer/docs edits don't bust the compile cache
COPY CMakeLists.txt ./
COPY src/ ./src/
COPY external/ ./external/
RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCYFTOOLS_WITH_TIFF=OFF \
        -DCYFTOOLS_WITH_CAIRO=OFF \
        -DCYFTOOLS_WITH_LDA=OFF \
        -DCYFTOOLS_WITH_OPENMP=OFF \
        -DCYFTOOLS_BUILD_TESTS=OFF \
    && cmake --build build -j "$(nproc)" --target cyftools

# ---------- stage 2: slim runtime ----------
# google-cloud-storage is only used when CYFVIEW_BUCKET is set (cloud big-file path);
# the import is lazy, so a purely-local container never touches it.
FROM ubuntu:24.04 AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip zlib1g ca-certificates \
    && pip3 install --break-system-packages --no-cache-dir google-cloud-storage \
    && rm -rf /var/lib/apt/lists/* /root/.cache \
    && useradd -m -u 10001 cyf
WORKDIR /app
COPY --from=build /src/build/cyftools /usr/local/bin/cyftools
COPY viewer/ /app/viewer/
ENV CYFTOOLS=/usr/local/bin/cyftools \
    CYFVIEW_CACHE=/tmp/cyfview-cache \
    CYFVIEW_MAX_CELLS=4000000 \
    PYTHONUNBUFFERED=1 \
    PORT=8080
EXPOSE 8080
USER cyf
CMD ["python3", "viewer/view.py"]
