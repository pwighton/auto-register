# Auto-Register Docker Image

This directory contains the Dockerfile for building the auto-register Docker image.

## Prerequisites

- Docker installed on your system
- A valid FreeSurfer license file

## Building the Image

1. Copy your FreeSurfer license file to this directory:
   ```bash
   cp /path/to/your/license.txt docker/.license
   ```

2. Build the Docker image:
   ```bash
   make
   ```

   This will create a Docker image tagged as `pwighton/areg`.

## Running the Container

### Interactive shell:
```bash
docker run -it --rm pwighton/areg
```

### Run auto_register.py:
```bash
docker run -it --rm \
  -v /path/to/data:/data \
  pwighton/areg \
  auto_register.py -s /data/output \
    --input-mode directory \
    --watch-directory /data/input \
    --disable-default-areg \
    --command /path/to/your/script.sh
```

## What's Included

- FreeSurfer 8.1.0
- Python 2.7 with required packages (numpy, nibabel, pydicom, etc.)
- dcm2niix for DICOM to NIfTI conversion
- mri_robust_register from FreeSurfer
- auto-register source code at `/opt/auto-register/src/`

## Notes

- The `.license` file is excluded from git via `.gitignore`
- The Python 2 environment is separate from any existing Python in the FreeSurfer container
- FreeSurfer binaries are available via the `$PATH` environment variable
