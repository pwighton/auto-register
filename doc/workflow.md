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

### Scanner host

#### Building the pulse sequences

TODO

#### Installing sequences on the host

TODO

### Registration laptop

#### Preparing the auto_register.py execution environment

TODO

#### Control room setup.

1. Connect the laptop to the scanner network via the switch behind the
host console.
1. Set the IP of the laptop to 192.168.TODO.TODO.

### Image registration machine

#### Building the ICE program

TODO

#### Installing the ICE program

TODO

## Execution ##

## Analysis ##
