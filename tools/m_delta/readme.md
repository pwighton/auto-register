# Compute Mdelta helper script

These helper scripts generate matrices that can be passed to auto-register with the `-prescrip` flag.

## `compute_MdeltaPrescription_v01.m`

Given a 'scout' image and:
- `pos`: A position (LPS) relative to the scout image, which you'd like to be the center of the new FOV
- `RO_dir`: A read-out direction (LPS) relative to the scout image, which you'd like to be the read-out direction of the new FOV
- `PE_dir`: A phase-encode direction (LPS) relative to the scout image, which you'd like to be the read-out direction of the new FOV
- `filestr`: name of the file to write the results to

Will produce:
- A text file containing a 4x4 matrix (in LPS) which can be passed to autoregister using the `-prescrip` flag.

Notes:
- Currently only works if the initial scout image was saggital at isocenter
- `compute_MdeltaPrescription_v01.m` was generously provided by Mukund Balasubramanian
- All position/direction vectors should be in LAS (Seimens and auto-register's native space)

## `compute_MdeltaPrescription_v01_RAS.m`

`compute_MdeltaPrescription_v01_RAS.m` is a wrapper around `compute_MdeltaPrescription_v01.m`.  It uses all the same parameters except all position/direction vectors (`pos`, `RO_dir`, `PE_dir`) should be in RAS rather than LPS.  This is a convinence script to make working with tools like freeview easier.

Notes:
- Currently only works if the initial scout image was saggital at isocenter
- All position/direction vectors should be in RAS

## `compute_MdeltaPrescription_fromSS_v01_RAS.m`

`compute_MdeltaPrescription_fromSS_v01_RAS.m` is a wrapper around `compute_MdeltaPrescription_v01_RAS.m`.

Given a 'scout' image and:
- `pos`: A position (RAS) relative to the scout image, which you'd like to be the center of the new FOV
- `SS_dir`: A slice selection direction (RAS) relative to the scout image which you'd like to be the slice select direction of the new FOV

Will produce:
- A text file containing a 4x4 matrix (in LPS) which can be passed to autoregister using the `-prescrip` flag.
- Given a desired `SS_dir`, there are an infinite number of `RO_dir`/`PE_dir` pairs that would satisfy the orthogonality constraints.  This script arbitrarily picks one.

## `Mdelta_RAS.bash`

`Mdelta_RAS.bash` is a convinence wrapper around `compute_MdeltaPrescription_fromSS_v01_RAS.m` which simplifies it's invocation

Example:
```
./Mdelta_RAS.bash "[1,2,3]" "[1,0,0]" mat01.txt
```

It can be copied out of this directory and placed in your working directory, simply change the `OCTAVE_SCRIPT_DIR` to point to the location of the `compute_MdeltaPrescription*.m` files

