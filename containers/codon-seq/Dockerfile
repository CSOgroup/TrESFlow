# syntax=docker/dockerfile:1.7

FROM --platform=linux/amd64 docker.io/library/python:3.12.13-slim-bookworm@sha256:4766d8b510c428e595d74b9cc5bbb2fae8e26316fffb4adc89908d79aacd58a2

ARG TARGETARCH
RUN test "${TARGETARCH}" = "amd64"

ADD --checksum=sha256:9b93e83b2aeb0dfb00d673a4482e9e22500469b92233438a2bbf188e0d08aa37 \
    https://github.com/exaloop/codon/releases/download/v0.16.3/codon-linux-x86_64.tar.gz \
    /tmp/codon-linux-x86_64.tar.gz
ADD --checksum=sha256:5684a5261a24fec7d3112349a1f0c00fbdba8237858aee2fcc6aa7beb4418e6f \
    https://github.com/exaloop/seq/releases/download/v0.11.3/seq-linux-x86_64.tar.gz \
    /tmp/seq-linux-x86_64.tar.gz
ADD --checksum=sha256:54a51198fca4e6eae44371277c525d3b7371cc969cd97e9f10ea6efc13acaf1a \
    https://deb.debian.org/debian/pool/main/p/pigz/pigz_2.6-1_amd64.deb \
    /tmp/pigz_2.6-1_amd64.deb
ADD --checksum=sha256:e82ba5d01929eafb8b9954606a3c38b0332a987c2b3432388b4ee7365e54deae \
    https://deb.debian.org/debian/pool/main/p/procps/libproc2-0_4.0.2-3_amd64.deb \
    /tmp/libproc2-0_4.0.2-3_amd64.deb
ADD --checksum=sha256:d9d0e75779cb79af869181f17b93c5c263a2b89cac6a0193c436160a4483ddc1 \
    https://deb.debian.org/debian/pool/main/p/procps/procps_4.0.2-3_amd64.deb \
    /tmp/procps_4.0.2-3_amd64.deb

RUN mkdir -p /opt/codon /opt/codon/lib/codon/plugins/seq \
    && tar -xzf /tmp/codon-linux-x86_64.tar.gz -C /opt/codon --strip-components=1 \
    && tar -xzf /tmp/seq-linux-x86_64.tar.gz -C /opt/codon/lib/codon/plugins/seq --strip-components=1
RUN mkdir -p /tmp/pigz-root /tmp/libproc-root /tmp/procps-root \
    && /usr/bin/dpkg-deb -x /tmp/pigz_2.6-1_amd64.deb /tmp/pigz-root \
    && /usr/bin/dpkg-deb -x /tmp/libproc2-0_4.0.2-3_amd64.deb /tmp/libproc-root \
    && /usr/bin/dpkg-deb -x /tmp/procps_4.0.2-3_amd64.deb /tmp/procps-root \
    && install -m 0755 /tmp/pigz-root/usr/bin/pigz /usr/local/bin/pigz \
    && ln -s pigz /usr/local/bin/unpigz \
    && install -m 0644 /tmp/libproc-root/lib/x86_64-linux-gnu/libproc2.so.0.0.1 /usr/lib/x86_64-linux-gnu/libproc2.so.0.0.1 \
    && ln -s libproc2.so.0.0.1 /usr/lib/x86_64-linux-gnu/libproc2.so.0 \
    && install -m 0755 /tmp/procps-root/bin/ps /usr/bin/ps \
    && /sbin/ldconfig \
    && rm -rf /tmp/pigz-root /tmp/libproc-root /tmp/procps-root \
    && rm -f /tmp/codon-linux-x86_64.tar.gz /tmp/seq-linux-x86_64.tar.gz \
        /tmp/pigz_2.6-1_amd64.deb /tmp/libproc2-0_4.0.2-3_amd64.deb \
        /tmp/procps_4.0.2-3_amd64.deb

ENV PATH="/opt/codon/bin:${PATH}"

COPY containers/codon-seq/smoke.codon /tmp/seq-smoke.codon
RUN test "$(codon --version)" = "0.16.3" \
    && grep -Eq '^version = "0\.11\.3"$' /opt/codon/lib/codon/plugins/seq/plugin.toml \
    && test "$(python3 --version)" = "Python 3.12.13" \
    && test "$(pigz --version)" = "pigz 2.6" \
    && ps -o pid= -p 1 >/dev/null \
    && codon run -plugin seq -release /tmp/seq-smoke.codon \
    && rm -f /tmp/seq-smoke.codon

LABEL org.opencontainers.image.title="TrESFlow Codon/Seq task runtime" \
      org.opencontainers.image.source="https://github.com/CSOgroup/TrESFlow" \
      org.opencontainers.image.version="codon-0.16.3_seq-0.11.3_python-3.12.13_pigz-2.6_procps-4.0.2"
