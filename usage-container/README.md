# Container workflow for HADDOCK3

This repository provides everything you need to build, run, and extend an **Apptainer/Singularity** container for HADDOCK3. [**HADDOCK3**](https://www.biorxiv.org/content/10.1101/2025.04.30.651432v1)(High Ambiguity Driven protein–protein Docking) is a flexible, information-driven software suite for modeling biomolecular complexes using experimental and theoretical restraints. **Docker**, **Singularity**, and **Apptainer** are containerization platforms that package applications and all their dependencies into lightweight, portable images, ensuring reproducible execution across different environments.

A ready-to-use Docker image for HADDOCK3 is published on the GitHub container registry, simply pull:
```bash
docker pull ghcr.io/haddocking/haddock3:2025.5.0
```
See the [HADDOCK3 container package](https://github.com/haddocking/haddock3/pkgs/container/haddock3) for details and version tags. This image serves as the canonical, versioned distribution of HADDOCK3 and can be used across development, CI/CD, and cloud environments.

To use this as the foundation for HPC-friendly SIF images, build **Apptainer** or **Singularity** containers directly from the Docker image in a single step. For example, to create an image:

```bash
# Build Apptainer image
apptainer build haddock3_mpi.sif docker://ghcr.io/haddocking/haddock3:2025.5.0

# Or, build with Singularity
singularity build haddock3_mpi.sif docker://ghcr.io/haddocking/haddock3:2025.5.0
  ```
**Note:** Refer to the [**usage.md**](https://github.com/Comp-era/HADDOCK3-Container/blob/main/docs/usage.md) and the example SLURM script available in the `scripts` folder for detailed instructions on how to run HADDOCK3 jobs in an HPC environment.
<pre> <strong>Version updates:</strong> The HADDOCK3 image is regularly updated. Please check the tags at <a href="https://github.com/haddocking/haddock3/pkgs/container/haddock3">ghcr.io/haddocking/haddock3</a> for the latest version. </pre>

---

### Build Your Own HADDOCK3 Container
Containerization techniques enable highly reproducible, customizable, and scalable scientific workflows. If you want to understand how the container is built under the hood from scratch or customize it for your own workflows, follow these detailed steps:
1. **Clone**

    ```bash
    git clone https://github.com/Comp-era/HADDOCK3-Container.git
    cd HADDOCK3-Container/recipe
    ```

    **Repository Structure**
    
    ```plaintext
    ├── apptainer_recipe/                  
    │   └── HADDOCK3.def              # Definition file
    ├── docs/                         # Documentation 
    │   ├── usage.md                  # Usage guide
    ├── scripts/                      # Scripts
    │   ├── slurm_run.sh              # Multi-node MPI run script
    ├── LICENSE                       # MIT License
    ├── README.md                     # Overview
    └── CONTRIBUTING.md              # Contribution guidelines
    ```

2. **Build** 

A definition file is a blueprint that tells the containerization platform how to build the container. It specifies the base OS, software packages, environment variables, and custom setup steps, ensuring your container is reproducible and tailored to your workflow.
A ready-to-use `HADDOCK3.def` is provided in the ([recipe/](https://github.com/Comp-era/HADDOCK3-Container/tree/main/recipe)) directory.

   ```bash
   # Apptainer
   apptainer build haddock3_mpi.sif haddock3_mpi.def

   # Singularity
   singularity build haddock3_mpi.sif haddock3_mpi.def
   ```

    **Tip:** To create your own definition files to layer in additional packages, alternative MPI variants, Python libraries, and any domain-specific utilities.Just simply modify the `%post` section of `HADDOCK3.def` before building.
    
3. **Verify&Run**

   - **Shell**: interactive access
     ```bash
     apptainer shell haddock3_cpu-mpi.sif
     ```
   - **HADDOCK3**: verify installation
     ```bash
     apptainer shell haddock3_cpu-mpi.sif haddock3 --version
     ```


4. **Download a pre-built image**

   You can also pull the pre-built apptainer image I created directly from GitHub’s container registry using the ORAS protocol:

   ```bash
   apptainer pull oras://ghcr.io/comp-era/haddock3:2025.06-v1.0-haddock3-mpi
   ```
---

##  Resources & Tutorials

- **Apptainer Installation & Usage**: Detailed installation instructions and usage examples can be found on the official docs: [apptainer.org/docs/admin/main/installation.html](https://apptainer.org/docs/admin/main/installation.html)
- **Official HADDOCK3 Tutorials**: Visit the Bonvin lab’s educational page for HADDOCK3 tutorials : [bonvinlab.org/education/HADDOCK3](https://www.bonvinlab.org/education/HADDOCK3/)
- **Source Code & Issues**: Explore the HADDOCK3 source code on GitHub: [github.com/haddocking/haddock3](https://github.com/haddocking/haddock3)

---

##  Requirements

- **Host**: Linux with Apptainer or Singularity installed (local machine or HPC environment)
- **Disk**: ≥ 2 GB free for building
- **Python**: 3.10+ (inside container)

---

##  Documentation

See the `docs/` folder:

- [**usage.md**](https://github.com/Comp-era/HADDOCK3-Container/blob/main/docs/usage.md) – Detailed usage and instructions to run Apptainer.

---

##  Contributing

Kindly refer to [CONTRIBUTING.md](CONTRIBUTING.md) , for contributions.



---

##  License

MIT License © 2025 Shantanu Khatri.

---

*Get ready for seamless reproducible workflow for HADDOCK3 runs!*

