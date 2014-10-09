AutoRegister experiment workflow documentation
==========

## Overview

This document explains the components and workflow required to perform
experiments for Phase One of the AutoRegister grant. These experiments
are designed to determine whether brain tumor size progression is more
reliably measured when MRI scans of a patient are acquired using the
same slice presciption across visits.

## Components

### Scanner host

The host computer configures and dispatches the pulse sequences
executed on the scanner. There are two pulse sequence programs
associated with AutoRegister: a diagnostic sequence for tumor imaging
that has been modified to allow changing the slice prescription, and a
scout sequence (AutoRegisterScout) that acquires a volume used to
determine the slice prescription.

### Registration laptop

The registration laptop is external to the Siemens scanning system,
and should run a modern Linux distribution. During the experiment, a
python script (auto_register.py) runs on the laptop. The script
receives image volumes from the image reconstruction machine,
interfaces with FreeSurfer tools to perform image volume registration,
and transmits slice prescription information to the scanner host.

### Image reconstruction machine

The image reconstruction machine reconstructs image volumes from the
raw data acquired by the scanner. The reconstruction (ICE) program for
the AutoRegisterScout sends image volumes to the registration laptop.

## Setup

### Preliminary

1. Clone the repo git@gitlab:ohinds/auto_register
1. Clone the repo git@gitlab:pwighton/vsend

### Scanner host

#### Building the pulse sequences

* AutoRegisterScout

1. Copy the directory seq/AutoRegisterScout from the auto\_register
repository to ${IDEA_BASE}/n4/pkg/MrServers/MrImaging/seq/
1. In the Siemens IDEA sde, change to the AutoRegisterScout sequence:
cs AutoRegisterScout
1. Issue the 'make sequence' command
ms

* Diagnostic sequence

1. Copy the directory seq/AutoRegisterApply from the auto\_register
repository to ${IDEA_BASE}/n4/pkg/MrServers/MrImaging/seq/
1. Follow the directions at the top of
AutoRegisterApply/AutoRegisterApply.h to use the transform application
tools inside the diagnostic sequence.

#### Installing sequences on the host

TODO

### Registration laptop

#### Preparing the auto_register.py execution environment

TODO

#### Control room setup.

1. Connect the laptop to the scanner network via the switch behind the
host console.
1. Set the IP of the laptop to 192.168.TODO.TODO.

### Image reconstruction machine

#### Building and installing the ICE program

* Follow the instructions in the Readme.txt file in the vsend repository.

## Execution ##

## Analysis ##
