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
