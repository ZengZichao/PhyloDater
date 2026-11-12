# PhyloDater Docker Image
# Multi-stage build for optimized image size

# =============================================================================
# Stage 1: Builder
# =============================================================================
FROM python:3.10-slim as builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
COPY pyproject.toml .
RUN pip install --upgrade pip && \
    pip install build && \
    pip install biopython ete3 ete4 dendropy pyyaml scipy numpy psutil matplotlib pandas

# =============================================================================
# Stage 2: Runtime - Minimal
# =============================================================================
FROM python:3.10-slim as runtime

LABEL maintainer="Zichao Zeng <zengzichao@sjtu.edu.cn>"
LABEL description="PhyloDater: A multi-software parallel platform for phylogenetic molecular dating"
LABEL version="0.1.0"
LABEL org.opencontainers.image.title="PhyloDater"
LABEL org.opencontainers.image.description="Multi-software parallel platform for phylogenetic molecular dating"
LABEL org.opencontainers.image.version="0.1.0"
LABEL org.opencontainers.image.licenses="MIT"

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    libopenblas-base \
    libgomp1 \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Set working directory
WORKDIR /app

# Copy application code
COPY . /app/

# Install PhyloDater
RUN pip install -e .

# Create directories for data and output
RUN mkdir -p /data /output

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PHYLODATER_LOG_LEVEL=INFO

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD phylodater --help || exit 1

# Default command
CMD ["phylodater", "--help"]

# =============================================================================
# Stage 3: Complete (with all external tools)
# =============================================================================
FROM runtime as complete

LABEL maintainer="Zichao Zeng <zengzichao@sjtu.edu.cn>"
LABEL description="PhyloDater Complete: Including all external bioinformatics tools"

# =============================================================================
# READ THIS BEFORE TRUSTING THE "complete" TAG  (review item INSTALL-002)
# -----------------------------------------------------------------------------
# The stage name and the LABEL above state an intent, not the content of this
# image.  As of v0.1.0 the `complete` stage ships THREE of PhyloDater's seven
# dating engines; the other four must be added by hand.
#
#   SHIPPED HERE (built from upstream sources further down this stage):
#     1. mcmctree   PAML 4.10.8 release tarball (abacus-gene/paml)
#     2. lsd2       IQ-TREE 2.3.4 Linux tarball -- LSD2 is the `--date` mode of
#                   iqtree2, so no separate LSD2 build is involved
#     3. pathd8     PATHd8.zip from www2.math.su.se, compiled with gcc -O3
#
#   NOT SHIPPED -- add them in a derived image (routes and caveats:
#   docs/EXTERNAL_TOOLS.md §2-§5; the same table in prose is §7):
#     4. treepl     no longer published on bioconda/conda-forge/PyPI; needs a
#                   source build (git clone blackrim/treePL; configure; make) or
#                   the historic conda 2.6.3 build repaired with patchelf.  The
#                   cmake/make/libnlopt-dev/libnlopt0 packages installed below
#                   exist precisely to make that build possible.
#     5. r8s        `pip install git+https://github.com/iTaxoTools/pyr8s.git`
#                   (NPRS only, range calibrations advisory) or the r8s binary
#                   from SourceForge
#     6. wlogdate   `pip install wlogdate` (provides launch_wLogDate.py)
#     7. mdcat      `pip install git+https://github.com/uym2/MD-Cat.git`
#
#   CONSEQUENCE: `--method treepl|r8s|wlogdate|mdcat` inside this image fails at
#   validate_environment() with a missing-executable error, and `phylodater check`
#   reports those engines as unavailable.  A seven-method run such as the one in
#   docs/TEST_RESULTS.md §1.1 was NOT executed in this image, so "complete" must
#   not be quoted as "all seven engines".
#
#   No new engine builds are added here on purpose: each of the four missing
#   routes needs its own source build and provenance record before it can be
#   claimed as supported.
# =============================================================================

# Install additional build dependencies for external tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    make \
    libnlopt-dev \
    libnlopt0 \
    && rm -rf /var/lib/apt/lists/*

# 1. Install PAML (for MCMCTree)
RUN cd /tmp && \
    wget -q https://github.com/abacus-gene/paml/releases/download/4.10.8/paml-4.10.8-linux-X86_64.tgz && \
    tar -xzf paml-4.10.8-linux-X86_64.tgz && \
    cp paml-4.10.8/bin/* /usr/local/bin/ && \
    rm -rf paml-4.10.8* && \
    echo "PAML installed successfully"

# 2. Install IQ-TREE2 (for LSD2)
RUN cd /tmp && \
    wget -q https://github.com/iqtree/iqtree2/releases/download/v2.3.4/iqtree-2.3.4-Linux.tar.gz && \
    tar -xzf iqtree-2.3.4-Linux.tar.gz && \
    cp iqtree-2.3.4-Linux/bin/iqtree2 /usr/local/bin/ && \
    rm -rf iqtree-2.3.4* && \
    echo "IQ-TREE2 installed successfully"

# 3. Install PATHd8
RUN cd /tmp && \
    wget -q http://www2.math.su.se/PATHd8/PATHd8.zip && \
    unzip -q PATHd8.zip && \
    cd PATHd8 && \
    gcc -O3 -o PATHd8 PATHd8.c -lm && \
    cp PATHd8 /usr/local/bin/ && \
    cd /tmp && \
    rm -rf PATHd8* && \
    echo "PATHd8 installed successfully"

# Verify installations
RUN echo "=== Installed Tools ===" && \
    which mcmctree && \
    which iqtree2 && \
    which PATHd8 && \
    echo "======================="

# =============================================================================
# Stage 4: Tutorial (with example data and Jupyter)
# =============================================================================
FROM complete as tutorial

LABEL maintainer="Zichao Zeng <zengzichao@sjtu.edu.cn>"
LABEL description="PhyloDater Tutorial: Complete environment with Jupyter notebooks"

# Install Jupyter and tutorial dependencies
RUN pip install --no-cache-dir \
    jupyter \
    jupyterlab \
    matplotlib \
    seaborn \
    pandas

# Copy tutorial notebooks
COPY tutorials/ /app/tutorials/
COPY examples/ /app/examples/

# Expose Jupyter port
EXPOSE 8888

# Set up Jupyter configuration
RUN jupyter notebook --generate-config && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.allow_root = True" >> /root/.jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.open_browser = False" >> /root/.jupyter/jupyter_notebook_config.py

# Default command for tutorial mode
CMD ["jupyter", "lab", "--allow-root", "--ip=0.0.0.0", "--port=8888", "--no-browser"]
