FROM python:3.12-slim-bookworm

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps bash \
    && rm -rf /var/lib/apt/lists/*
