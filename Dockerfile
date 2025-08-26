# 1) Match the R version in your lockfile
FROM rocker/r-ver:4.3.3

ENV DEBIAN_FRONTEND=noninteractive

# 2) System libs (IRMA and R)
RUN apt-get update && apt-get install -y --no-install-recommends \
    perl unzip \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libharfbuzz-dev libfribidi-dev libgit2-dev \
    ca-certificates curl git make g++ \
 && rm -rf /var/lib/apt/lists/*

# 3) Install renv + force stable CRAN mirror
RUN R -q -e "install.packages('renv', repos='https://cloud.r-project.org')"
ENV RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org

# 4) Keep renv library OUTSIDE /app so bind-mounts don’t hide it
ENV RENV_PATHS_LIBRARY=/opt/renv/library
RUN mkdir -p /opt/renv/library && chmod -R 777 /opt/renv

# 5) Restore dependencies from lockfile
WORKDIR /opt/app
COPY renv.lock ./
RUN R -q -e "renv::consent(provided=TRUE); renv::init(bare=TRUE); renv::restore()"

# 6) Rscripts
COPY main.R util.R README.md ./
COPY test.sh ./test.sh
RUN chmod +x ./test.sh
COPY test/ ./test/

# --- IRMA v1.2 (ZIP bundle that includes IRMA, IRMA_RES, LABEL, LABEL_RES) ---
ARG IRMA_URL=https://github.com/CDCgov/irma/releases/download/v1.2.0/flu-amd-202408.zip

RUN curl -L "$IRMA_URL" -o /tmp/irma.zip \
 && mkdir -p /opt/irma.tmp \
 && unzip -q /tmp/irma.zip -d /opt/irma.tmp \
 # move the single top-level dir inside the zip to /opt/irma
 && mv "$(find /opt/irma.tmp -mindepth 1 -maxdepth 1 -type d | head -n1)" /opt/irma \
 && rm -rf /tmp/irma.zip /opt/irma.tmp \
 # sanity checks: both IRMA and LABEL plus their *_RES must exist
 && test -x /opt/irma/IRMA && test -d /opt/irma/IRMA_RES \
 && test -x /opt/irma/LABEL && test -d /opt/irma/LABEL_RES

# Make both IRMA and LABEL discoverable
ENV PATH="/opt/irma:${PATH}"

# Optional smoke tests (don’t fail build if they print help with nonzero exit)
RUN IRMA --version 2>&1 | head -n 20 || true
RUN LABEL -h 2>&1 | head -n 20 || true

# 8) Default entrypoint
ENTRYPOINT ["Rscript","main.R"]
