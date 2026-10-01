"""Main FastAPI application entrypoint for HADDOCK3 Cloud GPU Service."""

from contextlib import asynccontextmanager
from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from server.config import settings
from server.routers import docking, results, system


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifecycle events: initialize workspace directory on startup."""
    settings.WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)
    yield


app = FastAPI(
    title="HADDOCK3 GPU Cloud Service",
    description=(
        "High-performance biomolecular docking control plane deployed on Google Cloud Run "
        "and powered by on-demand NVIDIA A100 serverless GPUs via Modal."
    ),
    version="3.0.0-gpu",
    lifespan=lifespan,
    docs_url="/docs",
    redoc_url="/redoc",
)

# Enable CORS for web portals and third-party scientific client integrations
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Register API v1 endpoints
app.include_router(docking.router, prefix="/api/v1")
app.include_router(results.router, prefix="/api/v1")
app.include_router(system.router, prefix="/api/v1")


@app.get(
    "/",
    tags=["Root"],
    summary="Service index and documentation links",
)
async def root() -> dict[str, str]:
    """Root endpoint welcoming users and pointing to interactive API specifications."""
    return {
        "title": "HADDOCK3 GPU Docking Cloud API",
        "description": "Submit biomolecules for accelerated docking on NVIDIA A100 GPUs.",
        "docs_url": "/docs",
        "redoc_url": "/redoc",
        "health_check": "/api/v1/system/health",
        "gpu_info": "/api/v1/system/gpu",
        "version": "3.0.0-gpu",
    }


@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """Catch unhandled errors and format clean JSON responses."""
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            "error": "InternalServerError",
            "message": str(exc),
            "path": str(request.url),
        },
    )


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "server.main:app",
        host=settings.HOST,
        port=settings.PORT,
        reload=settings.DEBUG,
    )
