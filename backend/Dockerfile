FROM python:3.12-slim

WORKDIR /app

# Install system dependencies if needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY app/ ./app/
COPY build_kmer_ext.py .

# Build the optional native k-mer accelerator.
#
# This is genuinely optional: kmer_native.py falls back to an equivalent pure-Python
# counter and the API reports which path ran, so a failed build must NOT fail the
# image -- hence `|| echo`. What it must not do is fail SILENTLY, which is what was
# happening before: the deployed service had no native path at all and reported
# "python" forever while the documentation advertised a measured speedup.
RUN pip install --no-cache-dir cython setuptools \
    && (python build_kmer_ext.py build_ext --inplace \
        && echo "native k-mer accelerator: BUILT" \
        || echo "native k-mer accelerator: NOT BUILT - falling back to Python") \
    && pip uninstall -y cython setuptools 2>/dev/null || true

# Baked at build time so /api/health can report exactly which commit is running.
# Render passes its commit SHA via RENDER_GIT_COMMIT at runtime; this ARG covers
# every other host.
ARG GIT_COMMIT=""
ENV GIT_COMMIT=$GIT_COMMIT

# Expose port (Render sets PORT env variable)
EXPOSE 8000

# Run uvicorn server
# Shell form so ${PORT} is expanded at runtime. Render assigns the port dynamically;
# the exec form above hardcoded 8000 and ignored it, so the container could bind a
# port the platform was not routing to.
CMD uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}
