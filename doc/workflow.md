AutoRegister experiment workflow documentation
==========

## Overview

This document explains the components and workflow required to perform
experiments for Phase One of the AutoRegister grant. These experiments
are designed to determine whether brain tumor size progression is more
reliably measured when MRI scans of a patient are acquired using the
same slice presciption across visits.

## Components

### Scanner host (part of the Siemens system)

The host computer configures and dispatches the pulse sequences
executed on the scanner. There are two pulse sequence programs
associated with AutoRegister: a diagnostic sequence for tumor imaging
that has been modified to allow changing the slice prescription, and a
scout sequence (AutoRegisterScout) that acquires a volume used to
determine the slice prescription.

### Registration laptop (External to the Siemens system)

The registration laptop is external to the Siemens scanning system,
and should run a modern Linux distribution. During the experiment, a
python script (auto_register.py) runs on the laptop. The script
receives image volumes from the image reconstruction machine,
interfaces with FreeSurfer tools to perform image volume registration,
and transmits slice prescription information to the scanner host.

### Image reconstruction machine (part of the Siemens system)

The image reconstruction machine reconstructs image volumes from the
raw data acquired by the scanner. The reconstruction (ICE) program for
the AutoRegisterScout sends image volumes to the registration
laptop. The ICE program used is the VSend functor, as developed by
Paul Wighton.

## Setup

### Common

To obtain the code for the AutoRegister prohject, clone the git
repository:

`git clone git@gitlab.nmr.mgh.harvard.edu:ohinds/auto_register`

### Scanner host

#### Building the pulse sequences

##### AutoRegisterScout

* Copy the directory `auto_register/seq` (from the git clone above)
  into the Siemens IDEA sequence development directory
  `n4\pkg\MrServers\MrImaging\seq`
* Start the `multi_idea` environment
* Start the SDE
* Change to the AutoRegisterScout sequence using `cs`, then choosing
  AutoRegisterScout
* Build the sequence using the `ms` command

##### Diagnostic sequence

TODO

### Preparing the auto_register.py execution environment on the registration laptop

* Ensure that `python 2.7` and `virtualenv` are installed.
* From the `auto_register` directory (from the git clone above),
execute:
`virtualenv arenv`
* Enter the virtual environment:
`source arenv/bin/activate`
* Install nibabel in the virtual environment:
`pip install git+git://github.com/nipy/nibabel`

### ICE program to send scout images

* Clone the VSend project: `git clone
  git@gitlab.nmr.mgh.harvard.edu:pwighton/vsend`
* Copy the vsend directory into the ICE program directory
  `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors`
* Rename the directory to `IceVSend`
* Start the `multi_idea` environment
* Change to the `IceVSend` directory
* Compile the ICE program by using `mi`

## Execution ##

### Sequence installation

#### AutoRegisterScout

* On the windows machine containing the IDEA environment, copy the files
  `n4\x86\prod\bin\AutoRegisterScout.dll` and
  `n4\i86\prob\lib\AutoRegisterScout.i86` to a USB drive, then to
  `C:\Medcom\MriCustomer\seq` on the scanner host.
* Add the AutoRegisterScout sequence to your protocol in the Protocol
  Editor tool on the host.

#### Diagnostic sequence

TODO

### Registration laptop

* Connect the laptop to the scanner network via the switch behind the
host console. If there are two switches, make sure you are connected
to the switch that is on the internal network with the host and image
reconstruction machines.
* Set the IP of the laptop to 192.168.2.5.

### ICE program installation

* On the windows machine containing the IDEA environment, copy the
  following files to a USB drive:
  * `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors\IceVSend\IceProgramVSend.evp`
  * `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors\IceVSend\xlinux\libIceVSend.so`
  * `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors\IceVSend\x86\IceVSend.dll`
  * `n4\pkg\MrServers\MrVista\Ice\IceIdeaFunctors\IceVSend\x86\IceVSend.evp`
* On the host, from the USB drive location:
`cp IceProgramVSend.evp C:\MedCom\MriCustomer\IceConfigurators
cp x86\IceVSend.dll C:\MedCom\bin
cp x86\IceVSend.evp C:\MedCom\bin
cp xlinux\libIceVSend.so C:\MedCom\MCIR\Med\lib`

### Scanning (sequence test scans only)

* On the registration laptop change to the `auto_register` directory,
  then enter the python environment by `source arenv/bin/activate`
* Execute the registration daemon:
`./auto_register.py -s <PATH_TO_DESIRED_IMAGE_DIRECTORY> -f`
* Run a localizer
* Run the AutoRegisterScout sequence (note that the daemon acknoledges
  receipt of the scout sequence)
* Move the phantom or subject to a different location
* Run the AutoRegisterScout sequence again (note that registration
  occurs and 'Transform ready to send' is printed by the daemon)
* Run the AutoRegisterScout a third time, noting that the image is
  aligned with the first scout scan, even though the positioning is
  the same as the second scan.

## Analysis

TODO
