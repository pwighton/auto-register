AutoRegister User Guide
==========

## Contents
1. Overview
1. Setup
1. Execution

## Overview

This document explains the components and workflow required to perform
experiments for Phase One of the AutoRegister grant. These experiments
are designed to determine whether brain tumor size progression is more
reliably measured when MRI scans of a patient are acquired using the
same slice prescription across visits.

## Components

### Scanner host (part of the Siemens system)

The host computer configures and dispatches the pulse sequences
executed on the scanner. AutoRegister requires a scout sequence
(AutoRegisterScout) that acquires a reference volume used to determine
the slice prescription. This sequence must be installed on the scanner
prior to conducting AutoRegister-enabled experiments.

### Registration laptop (External to the Siemens system)

The registration laptop is external to the Siemens scanning system,
and should run a modern Linux distribution. During the experiment, a
python script (auto_register.py) runs on the laptop. The script
receives the AutoRegisterScout volume(s) from the Siemens image
reconstruction machine, interfaces with FreeSurfer tools to perform
image volume registration, and transmits an affine transformation back
to the image reconstruction machine.

### Image reconstruction machine (part of the Siemens system)

The image reconstruction machine reconstructs image volumes from the
raw data acquired by the scanner. The reconstruction (ICE) program for
the AutoRegisterScout (called IceAutoRegisterInterface) sends image
volumes to the registration laptop and waits for the registration
laptop to send back an affine transformation. This transformation is
then passed along to the Siemens AutoAlign system, which allows it to
be applied on the scanner host to subsequent diagnostic scans.

## Setup

### Scanner host

#### Building the pulse sequences

##### AutoRegisterScout

* The scout has been tested for IDEA VB17 or VD13.
* Clone the pulse sequence:
  `git clone git@gitlab.nmr.mgh.harvard.edu:ohinds/auto_register_sequences`
* Copy the directory `auto_register_sequences` (from the git clone
  above) into the Siemens IDEA sequence development directory
  `n4\pkg\MrServers\MrImaging\seq`.
* Start the `multi_idea` environment
* Start the SDE
* Change to the AutoRegisterScout sequence using `cs`, then choosing
  AutoRegisterScout
* Build the sequence using the `ms` command
* Copy the files
  `n4\x86\prod\bin\AutoRegisterScout.dll` and
  `n4\i86\prob\lib\AutoRegisterScout.i86` to a USB drive, then to
  `C:\Medcom\MriCustomer\seq` on the scanner host.
* Add the AutoRegisterScout sequence to your protocol in the Protocol
  Editor tool on the host.
* Setup the AutoRegisterScout protocol with 2mm isotropic voxels, 128
  matrix size, and 2x GRAPPA acceleration. When TR and TE are
  minimized, this should yield an acquisition time of 1:56.???? CHECK THIS?????

### Preparing the auto_register.py execution environment on the registration laptop

* Clone the registration laptop software:
  `git clone git@gitlab.nmr.mgh.harvard.edu:ohinds/auto_register`
* Ensure that `python 2.7` and `virtualenv` are installed.
* From the `auto_register` directory (from the git clone above),
execute:
`virtualenv arenv`
* Enter the virtual environment:
`source arenv/bin/activate`
* Install nibabel in the virtual environment:
`pip install git+git://github.com/nipy/nibabel`

### ICE program

* Clone the VSend project: `git clone
  git@gitlab.nmr.mgh.harvard.edu:ohinds/ice_auto_register_interface`
* Copy the ice_auto_register_interface directory into the ICE program
  directory `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors`
* Rename the directory to `IceAutoRegisterInterface`
* Start the `multi_idea` environment
* Change to the `IceAutoRegisterInterface` directory
* Compile the ICE program by using `mi`
* Copy the following files to a USB drive:
  * `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors\IceAutoRegisterInterface\IceProgramAutoRegisterInterface.ipr`
  * `n4\linux\prod\lib\libIceAutoRegisterInterface.so`
  * `n4\x86\prod\bin\IceAutoRegisterInterface.dll`
  * `n4\x86\prod\bin\IceAutoRegisterInterface.evp`
* On the host, from the USB drive location, copy
  * IceProgramAutoRegisterInterface.ipr to `C:\Medcom\MriCustomer\ice\ohinds`
  * The rest of the files to `C:\Medcom\MriCustomer\ice`

## Execution ##

There are two types of scans the AutoRegister system is currently
useful for. The first is proof-of-concept scans, where two scan
sessions are simulated in a single session by removing the subject
from the bore and repositioning their head before a second set of
scans. The second is when there are actually multiple scan sessions,
where the AutoRegister scout scan needs to be saved in the first
session, and subsequent sessions require registration and application
of the resulting affine transformation. We'll refer to these different
types below as 'single session' and 'multi-session'.

### Registration laptop

* Connect the laptop to the scanner network via the switch behind the
host console. If there are two switches, make sure you are connected
to the switch that is on the internal network with the host and image
reconstruction machines.
* Set the IP of the laptop to 192.168.2.5. On ubuntu:
  * sudo service network-manager stop
  * sudo ifconfig eth0 192.168.2.5
* Change to the `auto_register` directory
* Enter the python environment by `source arenv/bin/activate`
* Change to the `auto_register` directory
* Execute the registration server:
  * For a single session scan execute the server in a mode where it
    takes the first received scout as the reference volume:
    `./auto_register.py -s <PATH_TO_DESIRED_IMAGE_DIRECTORY> -f`
  * For the first scan of a multi-session scan, also execute the
    server to store the first scout as the reference:
    `./auto_register.py -s <PATH_TO_DESIRED_IMAGE_DIRECTORY> -f`
  * For subsequent scans of a multi-session scan, execute the server
    and specify the location of the reference volume, which has been
    stored in the first session:
    `./auto_register.py -s <PATH_TO_DESIRED_IMAGE_DIRECTORY> -r <PATH_TO_REFERENCE_VOLUME>`

### Scanner host

#### First session

* Directly after running the diagnostic scan (the scan in which the
  tumor measurement will be conducted), add and open an
  AutoRegisterScout scan.
* Copy the slices of the diagnostic scan to the AutoRegisterScout. Do
  this by right clicking on the diagnostic scan while the
  AutoRegisterScout is open, chose 'Copy parameters', then choose
  'Slices'.
* Run the AutoRegisterScout sequence (note that the daemon acknowledges
  receipt of the scout sequence)

#### Subsequent sessions

* Directly before running the diagnostic scan, add and open an
  AutoRegisterScout scan.
* After the scout is acquired, the registration laptop will report
  that the scan was received, and register this scout volume to the
  reference volume.
* Once this is done, it will communicate the affine transformation
  back to the image reconstruction machine.
  * If this is successful, a cardinal symbol (a square cross made of a
    vertical and a horizontal line, with arrows on each end) will
    appear just to the left of the AutoRegisterScout sequence that was
    just run.
  * If this fails, a squiggly down arrow will appear farther to the
    left of the AutoRegisterScout sequence. At the moment, this
    happens periodically for unknown reasons. If this happens, just
    rerun the AutoRegisterScout sequence. Debugging this problem is
    still TODO.
* Add and open the diagnostic sequence.
* Make sure the AutoAlign region is set to 'Brain Basis'.
* Run the diagnostic sequence.

<!--  LocalWords:  AutoRegister AutoRegisterScout
 -->
