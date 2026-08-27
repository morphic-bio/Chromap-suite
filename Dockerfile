# syntax=docker/dockerfile:1.7

# Ubuntu 22.04, linux/amd64 manifest (glibc 2.35), resolved 2026-08-27.
# Keep the manifest digest immutable; update it deliberately with the release.
ARG UBUNTU_IMAGE="ubuntu:22.04@sha256:79676deb51ebb02885b0b9d33788e78a37cf1045ad79d1bb04c6a222c3556b3d"
ARG CHROMAP_SUITE_VERSION="1.0.1"
ARG CHROMAP_SUITE_REVISION="98a4da086f81b7cb159d8fe44efff2fb168e0785"

FROM ${UBUNTU_IMAGE} AS builder

ARG CHROMAP_SUITE_VERSION
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        libbz2-dev \
        libcurl4-gnutls-dev \
        libdeflate-dev \
        libhts-dev \
        liblzma-dev \
        libssl-dev \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build/chromap-suite
COPY . .

# Build every runtime entrypoint needed by standalone ATAC/ChIP/Hi-C modules
# and by the STAR Suite libchromap composition boundary.
RUN make clean \
    && make -j"$(nproc)" all chromap_lib_runner \
    && test "$(./chromap --version 2>&1)" = "${CHROMAP_SUITE_VERSION}" \
    && test "$(./chromap --upstream-version 2>&1)" = "0.3.3-r519" \
    && for binary in \
         chromap \
         rapidmacs \
         chromap_callpeaks \
         chromap_lib_runner \
         chromap_atac_spill_materializer; do \
         test -x "${binary}"; \
       done

RUN install -d \
        /opt/chromap-suite/bin \
        /opt/chromap-suite/include/chromap-suite \
        /opt/chromap-suite/include/htslib \
        /opt/chromap-suite/include/rapidmacs \
        /opt/chromap-suite/lib \
    && install -m 0755 \
        chromap \
        rapidmacs \
        chromap_callpeaks \
        chromap_lib_runner \
        chromap_atac_spill_materializer \
        /opt/chromap-suite/bin/ \
    && install -m 0644 \
        libchromap.a \
        third_party/rapidmacs/lib/librapidmacs.a \
        /opt/chromap-suite/lib/ \
    && find src -maxdepth 1 -type f \( -name '*.h' -o -name '*.hpp' \) \
        -exec install -m 0644 '{}' /opt/chromap-suite/include/chromap-suite/ ';' \
    && cp -a third_party/htslib/htslib/. /opt/chromap-suite/include/htslib/ \
    && cp -a third_party/rapidmacs/include/rapidmacs/. \
        /opt/chromap-suite/include/rapidmacs/

# Prove that the installed prefix, rather than the source-tree include layout,
# is sufficient for a downstream libchromap consumer to compile and link.
RUN g++ -std=c++11 -fopenmp \
        -I/opt/chromap-suite/include/chromap-suite \
        -I/opt/chromap-suite/include \
        tests/container_sdk_link_smoke.cc \
        /opt/chromap-suite/lib/libchromap.a \
        /opt/chromap-suite/lib/librapidmacs.a \
        -lhts -lm -lz -lpthread -ldl -lcurl -lcrypto -lbz2 -llzma -ldeflate \
        -o /tmp/chromap_sdk_link_smoke \
    && /tmp/chromap_sdk_link_smoke

FROM ${UBUNTU_IMAGE} AS runtime

ARG CHROMAP_SUITE_VERSION
ARG CHROMAP_SUITE_REVISION
ENV DEBIAN_FRONTEND=noninteractive

LABEL org.opencontainers.image.title="Chromap Suite" \
      org.opencontainers.image.description="Chromap, libchromap, and RapidMACS for bulk and single-cell chromatin workflows" \
      org.opencontainers.image.version="${CHROMAP_SUITE_VERSION}" \
      org.opencontainers.image.revision="${CHROMAP_SUITE_REVISION}" \
      org.opencontainers.image.source="https://github.com/morphic-bio/Chromap-suite" \
      org.opencontainers.image.url="https://github.com/morphic-bio/Chromap-suite/releases/tag/v${CHROMAP_SUITE_VERSION}" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.base.name="ubuntu:22.04" \
      org.opencontainers.image.base.digest="sha256:79676deb51ebb02885b0b9d33788e78a37cf1045ad79d1bb04c6a222c3556b3d"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libbz2-1.0 \
        libcurl4 \
        libdeflate0 \
        libgomp1 \
        libhts3 \
        liblzma5 \
        libssl3 \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/chromap-suite /opt/chromap-suite

ENV PATH="/opt/chromap-suite/bin:${PATH}" \
    CHROMAP_TEMP_DIR="/tmp"

WORKDIR /data

RUN test "$(chromap --version 2>&1)" = "${CHROMAP_SUITE_VERSION}" \
    && for binary in \
         chromap \
         rapidmacs \
         chromap_callpeaks \
         chromap_lib_runner \
         chromap_atac_spill_materializer; do \
         command -v "${binary}"; \
         ! ldd "$(command -v "${binary}")" | grep -q "not found"; \
       done

CMD ["chromap", "--help"]
