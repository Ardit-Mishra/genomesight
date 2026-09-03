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

# Expose port (Render sets PORT env variable)
EXPOSE 8000

# Run uvicorn server
# Shell form so ${PORT} is expanded at runtime. Render assigns the port dynamically;
# the exec form above hardcoded 8000 and ignored it, so the container could bind a
# port the platform was not routing to.
CMD uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}
