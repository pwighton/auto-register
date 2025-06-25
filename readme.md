# AutoRegister

Custom slice prescription for Seimens MR scanners.  Send volumes to a remote computer, register them, and send the result back to the scanner and overwrite the [AutoAlign](https://www.siemens-healthineers.com/en-us/magnetic-resonance-imaging/options-and-upgrades/clinical-applications/autoalign) matrix.

Volumes can be sent to a remote computer using either [vsend](https://github.com/mriMotionCorrection/vsend) or the DICOM export option of the [Seimens BOLD dot addin](https://www.siemens-healthineers.com/en-us/magnetic-resonance-imaging/options-and-upgrades/clinical-applications/bold-3d-evaluation)  The [IceAASetter ICE functor](https://github.com/pwighton/IceAASetter/) queries the auto-regsiter server, retrives the latest registration matrix and overwrite the AutoAlign Matrix.

After the AutoAlign matrix has been set by the [IceAASetter ICE functor](https://github.com/pwighton/IceAASetter/), any sequence with AutoAlign set to `Head -> Brain` will have its FOV be transformed by the Auto Align matrix.

## Prereqs

- A computer that can be connected to the scanner with:
  - Ubuntu, or equivalent Linux distribution
  - Python 2.7 (conda recommended)
  - FreeSurfer.  Any version will do as long as it has
    - [mri_robust_register](https://surfer.nmr.mgh.harvard.edu/fswiki/mri_robust_register)
    - [mri_synthstrip](https://surfer.nmr.mgh.harvard.edu/docs/synthstrip/) if using the `--synthstrip` flag.
  - Samba if using the DICOM export option from the Seimens BOLD dot addin (`--input-mode directory`)
  - Git

## Installation

1) Clone this repo.

```
git clone git@github.com:pwighton/auto-register.git
```

2) Configure a python environment.

Conda is recommended, but any python enviornment manager should work.

```
conda create --name areg python=2.7
conda activate areg
pip install -r ./requirements.txt
```

3) Install the and configure [IceAASetter ICE functor](https://github.com/pwighton/IceAASetter/) on your Seimens scanner

4) Determine how data will be transfered to the auto-register computer

- If using vsend, Install and configure the [VSend ICE functor](https://github.com/pwighton/VSend) on your Seimens scanner
- If using the DICOM export option of the Seimens BOLD dot addin, [configure the BOLD dot addin](doc/bold-dotaddin-config-instructions.md) on your Seimens scanner

## Run

1) Connect the auto-register computer to the internal Seimens network and set the IP address accordingly.

```
sudo ifconfig eno2 192.168.2.5
```

- `192.168.2.5` is typically used, but this IP address must match:
  - The IP address vsend has been configured with; or
  - The IP address the BOLD dot addin has been configured with
- The interface (`eno2`) will vary.  Run `ifconfig` to find the propoper interface
- If you are using a connection manager, disable it so that it doesn't override the `ifconfig` command

2) Confirm the autoregister computer can ping the console computer and mars

```
ping 192.168.2.1
ping 192.168.2.2
```

3) Activate the FreeSurfer environment and the auto-register python environment

```
export FREESURFER_HOME /path/to/your/freesurfer/installation
source $FREESURER_HOME/SetUpFreeSurfer.sh
conda activate areg
```

4) Make a working diretory for auto-register for this session

```
mkdir ~/data/areg-dir-todays-date
```

5) Start auto-register

Run `src/auto_register.py --help` for a full list of options

### If using vsend

```
src/auto_register.py \
  -s ~/data/areg-dir-todays-date \
  -f
```

### If using the BOLD dot addin

```
src/auto_register.py \
  -s ~/data/areg-dir-todays-date \
  -f \
  --input-mode directory \
  --watch-directory /BOLD_DICOM
```

- `--watch-directory /BOLD_DICOM` should correspond to the directory the bold dot addin has been configured for (step #4 in the [BOLD dot addin configuration instructions](doc/bold-dotaddin-config-instructions.md))

The Bold dot addin transfer method assumes a 1-to-1 mapping between dicoms files and image volumes (i.e. 3d scans).  It's ok if a single sequence generates multiple volumes, as long as the 1-to-1 mapping between dicom files and volumes hold.  Use the `--dicom-filter` flag and point it to a json file to filter out volumes you'd like autoregsiter to ignore.  For example the, the json file [`dicom-filters/memprage-rms.json`](dicom-filters/memprage-rms.json) can be used to filter out all but the RMS volume of a multi-echo MPRAGE sequence.  Note: this filter assumes the sequence is named "MEMPRAGE_4mm_scout" on the scanner.  Edit accordingly.

## References

todo


