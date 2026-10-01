# HADDOCK3 GPU Cloud Service

A production-grade, asynchronous REST API enabling external clients to submit molecular structures (proteins, peptides, nucleic acids) for accelerated biomolecular docking powered by **NVIDIA A100 / T4 GPUs** via Modal, deployed on **Google Cloud Run**.

---

## 🏛️ Architecture Overview

The system uses a decoupled control plane / compute plane architecture:

1. **Control Plane (Google Cloud Run)**:
   - Serverless CPU container running FastAPI and Uvicorn.
   - Handles client HTTP requests, multipart file uploads, validation, and status tracking.
   - Automatically scales to zero when idle, minimizing idle hosting costs.
2. **Compute Plane (Modal Serverless Cloud GPUs)**:
   - On-demand NVIDIA A100 (or T4/A10G/H100) GPU workers.
   - Pre-baked CUDA 13.0 container with PyTorch, OpenMM GPU, and full HADDOCK3 suite.
   - Dispatched asynchronously via `.spawn()`, executing docking workflows (`topoaa`, `rigidbody`, `openmm`, `clustfcc`, `caprieval`).
3. **Artifact Persistence**:
   - Ranked docked structures and CAPRI evaluations are streamed back and exposed via REST endpoints for interactive inspection and download.

---

## 🚀 Quick Start (Local Development)

### 1. Requirements
Ensure Modal is authenticated:
```bash
modal auth
```

### 2. Run Server Locally
```bash
./server/run_server.sh
```
Or directly with Uvicorn:
```bash
uvicorn server.main:app --host 0.0.0.0 --port 8080 --reload
```

Interactive documentation is available at:
- Swagger UI: `http://localhost:8080/docs`
- ReDoc: `http://localhost:8080/redoc`

---

## 🧪 Testing

Run the full integration test suite:
```bash
pytest server/tests/test_server.py -v
```

---

## 📡 API Reference

### 1. Submit Docking Job
```bash
POST /api/v1/docking/submit
Content-Type: multipart/form-data
```
**Parameters**:
- **Molecule 1**:
  - `mol1_file`: PDB or mmCIF structure file (`.pdb`, `.ent`, `.cif`). Required.
  - `mol1_chain`: Which chain to use (default: `"All"`, or specific chain ID e.g. `"A"`).
  - `mol1_kind`: Molecule kind (default: `"Protein or Protein-Ligand"`, `"DNA"`, `"RNA"`, `"Small Molecule"`).
  - `mol1_coarse_grain`: Convert to Martini coarse-grained representation (`true`/`false`, default: `false`).
  - `mol1_cyclic_peptide`: Is it a cyclic peptide? (`true`/`false`, default: `false`).
- **Molecule 2**:
  - `mol2_file`: PDB or mmCIF structure file (`.pdb`, `.ent`, `.cif`). Required.
  - `mol2_chain`: Which chain to use (default: `"All"`, or specific chain ID e.g. `"B"`).
  - `mol2_kind`: Molecule kind (default: `"Protein or Protein-Ligand"`).
  - `mol2_coarse_grain`: Convert to Martini coarse-grained representation (`true`/`false`, default: `false`).
  - `mol2_cyclic_peptide`: Is it a cyclic peptide? (`true`/`false`, default: `false`).
- **Restraints & Compute Settings**:
  - `restraints`: Optional interaction restraint file (`.tbl` or `.act`).
  - `sampling`: Number of rigid-body models (default: `100`).
  - `refinement`: Number of flexible refinement models (default: `20`).
  - `gpu_type`: Target GPU accelerator (`A100`, `T4`, `A10G`, `H100`). Default: `A100`.
  - `job_name`: Optional title.

**Example `curl`**:
```bash
curl -X POST "http://localhost:8080/api/v1/docking/submit" \
  -F "mol1_file=@e2a.pdb" \
  -F "mol1_chain=A" \
  -F "mol1_cyclic_peptide=false" \
  -F "mol2_file=@hpr.pdb" \
  -F "mol2_chain=All" \
  -F "mol2_cyclic_peptide=false" \
  -F "restraints=@ti.tbl" \
  -F "sampling=100" \
  -F "refinement=20" \
  -F "gpu_type=A100"
```

**Response (`202 Accepted`)**:
```json
{
  "job_id": "8d3e24b1-8729-4d64-839e-21ef186c7d91",
  "status": "running",
  "gpu_type": "A100",
  "modal_call_id": "fc-01HF...",
  "status_url": "/api/v1/docking/8d3e24b1-8729-4d64-839e-21ef186c7d91/status",
  "results_url": "/api/v1/results/8d3e24b1-8729-4d64-839e-21ef186c7d91",
  "created_at": "2026-10-01T03:50:00+00:00"
}
```

---

### 2. Poll Job Status
```bash
GET /api/v1/docking/{job_id}/status
```
**Example `curl`**:
```bash
curl -X GET "http://localhost:8080/api/v1/docking/8d3e24b1-8729-4d64-839e-21ef186c7d91/status"
```

**Response**:
```json
{
  "job_id": "8d3e24b1-8729-4d64-839e-21ef186c7d91",
  "status": "running",
  "gpu_type": "A100",
  "current_stage": "05_clustfcc",
  "completed_stages": ["00_topoaa", "01_rigidbody", "02_seletop", "03_flexref", "04_emref"],
  "total_stages": 9,
  "elapsed_seconds": 45.2,
  "error_message": null
}
```

---

### 3. Retrieve Results
```bash
GET /api/v1/results/{job_id}
```
**Response**:
```json
{
  "job_id": "8d3e24b1-8729-4d64-839e-21ef186c7d91",
  "status": "completed",
  "gpu_type": "A100",
  "execution_seconds": 68.4,
  "clusters": [
    {
      "cluster_rank": 1,
      "cluster_id": 1,
      "model_count": 4,
      "haddock_score": -89.42,
      "fnat": 0.82,
      "irmsd": 1.25,
      "lrmsd": 2.10,
      "dockq": 0.81
    }
  ],
  "best_models": ["cluster1_1.pdb", "cluster1_2.pdb"],
  "download_urls": {
    "cluster1_1.pdb": "/api/v1/results/8d3e24b1-8729-4d64-839e-21ef186c7d91/download/cluster1_1.pdb"
  }
}
```

---

### 4. Download PDB Artifact
```bash
GET /api/v1/results/{job_id}/download/{filename}
```
Streams the docked PDB structure file directly.

---

### 5. Cancel Job
```bash
POST /api/v1/docking/{job_id}/cancel
```
Terminates the remote Modal cloud GPU execution.

---

## ☁️ Google Cloud Run Deployment

### 1. Build and Push Container Image to Google Artifact Registry
```bash
# Set your GCP Project ID
export PROJECT_ID="your-gcp-project-id"
export REGION="us-central1"
export IMAGE_TAG="$REGION-docker.pkg.dev/$PROJECT_ID/haddock3/gpu-server:latest"

# Build using Google Cloud Build
gcloud builds submit --tag "$IMAGE_TAG" -f server/Dockerfile .
```

### 2. Deploy to Google Cloud Run
Pass your Modal token credentials as environment variables or Google Secret Manager references:
```bash
gcloud run deploy haddock3-gpu-server \
  --image "$IMAGE_TAG" \
  --platform managed \
  --region "$REGION" \
  --allow-unauthenticated \
  --memory 2Gi \
  --cpu 2 \
  --min-instances 0 \
  --max-instances 10 \
  --set-env-vars "MODAL_TOKEN_ID=$MODAL_TOKEN_ID,MODAL_TOKEN_SECRET=$MODAL_TOKEN_SECRET"
```

Once deployed, Google Cloud Run returns a public HTTPS URL (e.g. `https://haddock3-gpu-server-xyz.a.run.app`) ready to serve external users worldwide!
