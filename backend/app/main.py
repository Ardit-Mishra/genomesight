from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.endpoints import router as api_router

app = FastAPI(
    title="GenomeSight API",
    description="High-performance bioinformatics REST API for sequence analysis, alignment, translation, and codon usage.",
    version="2.0.0"
)

# Configure CORS for frontend access (Vercel / Localhost)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for seamless developer experience & production flexibility
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(api_router)

@app.get("/")
def root():
    return {
        "project": "GenomeSight API",
        "status": "online",
        "docs_url": "/docs",
        "health_endpoint": "/api/health"
    }
