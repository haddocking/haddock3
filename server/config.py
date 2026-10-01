"""Configuration settings for the HADDOCK3 FastAPI Cloud Service."""

import os
from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings with environment variable overrides."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # Server configuration (Cloud Run provides PORT automatically)
    HOST: str = "0.0.0.0"
    PORT: int = int(os.environ.get("PORT", "8080"))
    DEBUG: bool = False

    # Modal cloud compute configuration
    MODAL_APP_NAME: str = "haddock3-gpu-suite"
    MODAL_ENVIRONMENT: str = "main"
    DEFAULT_GPU_TYPE: str = "A100"  # Options: T4, A10G, A100, H100

    # Local workspace for staging uploads before dispatching
    BASE_DIR: Path = Path(__file__).resolve().parent
    WORKSPACE_DIR: Path = BASE_DIR / "workspace"
    MAX_UPLOAD_SIZE_MB: int = 50

    # CORS configuration
    ALLOWED_ORIGINS: list[str] = ["*"]


settings = Settings()
