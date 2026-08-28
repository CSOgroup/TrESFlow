# Codon/Seq task runtime

This build is deliberately limited to `linux/amd64`. Codon 0.16.3 and Seq
0.11.3 do not publish Linux ARM64 release artifacts.

Build and smoke-test the local image from the repository root:

```bash
docker build --platform linux/amd64 \
  --file containers/codon-seq/Dockerfile \
  --tag tresflow-codon-seq:phase2-local .
docker run --rm --platform linux/amd64 \
  tresflow-codon-seq:phase2-local codon --version
```

Every remote build input and the base image are pinned in both the Dockerfile
and `runtime-manifest.json`. The Docker build also executes a real Seq
`import bio` smoke test.

This image must not be published yet. Codon 0.16.3's BSL 1.1 change date was
2026-05-01 and its change license is Apache-2.0. The upstream Seq v0.11.3 tag,
however, contains no license file or explicit redistribution grant. The HCC
Conda recipe labels its package Apache-2.0, but that downstream assertion is
not a substitute for an upstream grant. Production process container
directives therefore remain intentionally unset until Seq redistribution is
clarified and the built image is published under an immutable registry digest.
