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
    # The CORS spec forbids "*" together with credentials, and browsers reject the
    # pair outright. Starlette papers over it by echoing the request origin, which
    # means the previous config was effectively "any origin, WITH credentials" --
    # strictly more permissive than intended.
    #
    # This API is public and stateless: it reads no cookies and issues no session.
    # So the honest configuration is a genuine wildcard with credentials OFF, which
    # is both what the service needs and what browsers will actually honour.
    allow_origins=["*"],
    allow_credentials=False,
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
