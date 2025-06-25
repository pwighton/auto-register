# BOLD dot addin configuration instructions

1) On the auto-register computer, determine the user that will be used to run auto-register, and log in as that user

- This user needs to be able to run `sudo`
- In this example, the user will be called `paul`

2) Install samba

```
sudo apt-get install --upgrade samba
```

3) Run `smbpasswd`, set the samba password for the user

```
sudo smbpasswd -a paul
```

- In this example, the samba password for the user `paul` will be `paul123`
- Replace the user `paul` with the user from step #1
- There will be 3 password prompts
- The first password will be for sudo (enter the user's password)
- The second password will be the samba password to set for user `paul` (`paul123`)
- The third password will be to verify the second password (`paul123`)

4) Make the directory for incoming BOLD dicoms:

In this example, it will be called `/BOLD_DICOM`

```
sudo mkdir /BOLD_DICOM
sudo chown paul /BOLD_DICOM
sudo chmod -R 773 /BOLD_DICOM
```
- replace `paul` with the user from step #1

5) Configure samba 

- Edit the included `smb.conf` file and change the `path` from `/BOLD_DICOM` to the directory chosen in step #4
- Backup your current `smb.conf` file and replace it with the `smb.conf` in this bundle

```
sudo cp /etc/samba/smb.conf /etc/samba/smb.conf.backup.20250113
sudo cp /path/to/this/new/smb.conf /etc/samba/smb.conf
```

6) Connect the auto-register computer to the internal Siemens scanner network, and assign it an internal IP address

```
sudo ifconfig eno2 192.168.2.5
```

- replace `eno2` with the name of the ethernet interface (run `ifconfig` to see a list)
- confirm connection by pinging the console computer and mars
- `ping 192.168.2.1`
- `ping 192.168.2.2`

6) Start or Restart samba

```
sudo systemctl restart smbd nmbd
```

7) Confirm that samba is running

```
sudo systemctl status smbd nmbd
```

8) Confirm the `BOLD_DICOM` directory is being shared by samba

```
smbclient -L localhost -U paul
```

- replace `paul` with the user from step #1

should give something like:

```
	Sharename       Type      Comment
	---------       ----      -------
	BOLD_DICOM      Disk      Samba Share for BOLD dot-addin dicoms
	IPC$            IPC       IPC Service (Samba 4.15.13-Ubuntu)
```

9) Configure the BOLD dot-addin accordingly

On the scanner, add the BOLD dot-addin to the desired BOLD sequence, and configure:

- Target host: `\\192.168.2.5` (the IP used in step #6)
- Target Directory: \BOLD_DICOM (the directory used in step #4)
- User Name: WORKGROUP\paul (the user from step #1)
- Password: paul123 (the password from step #3)

10) Run a sequence with the BOLD dot-addin.  The DICOMs should appear in the `\BOLD_DICOM` folder of the auto-register computer once the scan completes.
