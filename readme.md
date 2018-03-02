# AutoRegister

## Prereqs

- Find a computer to run AutoRegister
- Install ubuntu
- Install git
- Install docker
- Install make
- Clone this repo

## Setup

- Build the areg docker container
  - `make docker`
- Install the AutoRegsiter sequence and Ice functor
  - Latest version of the binaries can be found here:
    - `/autofs/cluster/gerenuk/pwighton/AutoRegP2/binaries/`
- Connect the autoregsiter computer to the internal MR network (with IP 192.168.2.5)
  - Disconnect wifi and eth interfaces (from the ubuntu GUI)
  - Determine the network interface you are using
    - Traditionally, this has been `eth0` on ubuntu
    - Now it's `enx9cebe8347c61` on my machine for some reason
    - More info [here](https://www.freedesktop.org/wiki/Software/systemd/PredictableNetworkInterfaceNames/)
  - Open a terminal and run
    - `sudo ifconfig enx9cebe8347c61 192.168.2.5`
      - (replace `enx9cebe8347c61` with the network interface you are using)
  - Confirm the connection by trying to ping the MR console machine (and/or mars)
    - `ping 192.168.2.1`
    - `ping 192.168.2.2`
- Launch areg
  - Areg assumes imaging data is sent to areg on port `15000`
  - Areg assumes transform data is sent back to the scanner on port `15001`
  - Determine the location of this repo
    - e.g `/home/paul/lcn/git/areg/auto_register`
  - If this is a single session (or first in a series)
    - create a directory to store the session data
      - `mkdir -p /tmp/areg-data/test-subject`
    - launch areg, pointing to this dir
      - `docker run -it --rm \
           -v /home/paul/lcn/git/areg/auto_register:/areg \
           -v /tmp/areg-data:/tmp/areg-data \
           -w /areg \
           -u ${UID} \
           -p 15000:15000 \
           -p 15001:15001 \
           areg /areg/src/auto_register.py -s /tmp/areg-data/test-subject -f -H 127.0.0.1`
  - If this is a repeat session
    -

## Oliver's README (old)
Automatically register subsequent acquisitions to a reference.

For installation and usage instructions, see doc/workflow.md


TODO (oliver):

- Debug sometimes AA matrix doesn't take
  - Pipe service error. have logs
- Tests (nose)
