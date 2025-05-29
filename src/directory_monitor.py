#!/usr/bin/env python2

"""Directory monitoring functionality for ImageReceiver.
Monitors a directory for new DICOM files and converts them to NIfTI.
"""

import os
import glob
import time
import subprocess
import threading
import tempfile
import shutil
from collections import defaultdict

class DirectoryMonitor(object):
    """Monitor a directory for new DICOM files, convert to NIfTI, and make
    filenames available through the same interface as TCP image reception.
    """

    def __init__(self, args):
        """Initialize directory monitoring with command line arguments."""
        self.watch_directory = args.watch_directory
        self.save_location = args.save_directory  # Uses existing -s argument
        self.poll_interval = getattr(args, 'poll_interval', 2.0)
        self.stabilize_wait = getattr(args, 'stabilize_wait', 3.0)
        
        # Thread control
        self.mutex = threading.Lock()
        self.filename_stack = []
        self.should_stop = False
        self.monitor_thread = None
        
        # Track files we've seen to avoid reprocessing
        self.known_files = set()
        self.file_sizes = {}  # Track file sizes for stability detection
        self.processing_files = set()  # Files currently being processed
        
        # Counter for output filenames
        self.save_volume_index = 0
        
        # Initialize known files (ignore existing files at startup)
        self._scan_existing_files()
    
    def check_environment(self):
        """Make sure that our environment is able to execute dcm2niix
        """
        cmd = ['dcm2niix']
        try:
            check_proc = subprocess.Popen(cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
            out, err = check_proc.communicate()
            if (check_proc.returncode == 0 ):
                return True
            else:
                print "Unexpected output while executing %s" % cmd
                return False
        except(OSError):
            print "Can't find %s, make sure it's in the $PATH?" % cmd
            return False

    def _scan_existing_files(self):
        """Scan directory for existing DICOM files to ignore at startup."""
        print "Scanning existing DICOM files to ignore..."
        existing_files = set()
        
        # Use os.walk for Python 2 compatibility (no recursive glob)
        for root, dirs, files in os.walk(self.watch_directory):
            for filename in files:
                if filename.lower().endswith(('.dcm', '.ima')):
                    filepath = os.path.join(root, filename)
                    existing_files.add(filepath)
        
        self.known_files.update(existing_files)
        print "Found %d existing DICOM files to ignore" % len(existing_files)
    
    def start(self):
        """Start the directory monitoring thread."""
        if self.monitor_thread is not None:
            raise RuntimeError('Directory monitor already running')
        
        self.should_stop = False
        self.monitor_thread = threading.Thread(target=self._monitor_loop)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        print "Directory monitor started, watching: %s" % self.watch_directory
    
    def stop(self):
        """Stop the directory monitoring thread."""
        if self.monitor_thread is not None:
            self.should_stop = True
            self.monitor_thread.join()
            self.monitor_thread = None
        print "Directory monitor stopped"
    
    def get_next_filename(self):
        """Return the next available NIfTI filename, or None if no new files."""
        filename = None
        self.mutex.acquire()
        if len(self.filename_stack) > 0:
            filename = self.filename_stack.pop()
        self.mutex.release()
        return filename
    
    def is_running(self):
        """Check if the monitor is currently running."""
        return self.monitor_thread is not None and not self.should_stop
    
    def _monitor_loop(self):
        """Main monitoring loop - runs in separate thread."""
        print "Directory monitoring loop started"
        
        while not self.should_stop:
            try:
                self._scan_for_new_files()
                time.sleep(self.poll_interval)
            except Exception as e:
                print "Error in directory monitoring: %s" % str(e)
                # Continue monitoring despite errors
                time.sleep(self.poll_interval)
    
    def _scan_for_new_files(self):
        """Scan directory for new DICOM files."""
        current_files = set()
        
        # Use os.walk for Python 2 compatibility (no recursive glob)
        for root, dirs, files in os.walk(self.watch_directory):
            for filename in files:
                if filename.lower().endswith(('.dcm', '.ima')):
                    filepath = os.path.join(root, filename)
                    current_files.add(filepath)
        
        # Find truly new files (not in known_files and not being processed)
        new_files = current_files - self.known_files - self.processing_files
        
        for filepath in new_files:
            if self._is_file_stable(filepath):
                # File is stable, start processing
                self.processing_files.add(filepath)
                self.known_files.add(filepath)
                
                # Process in separate thread to avoid blocking monitoring
                process_thread = threading.Thread(
                    target=self._process_dicom_file, 
                    args=(filepath,)
                )
                process_thread.daemon = True
                process_thread.start()

    def _is_file_stable(self, filepath):
        """Check if file size has stabilized (file write is complete)."""
        try:
            current_size = os.path.getsize(filepath)
        except OSError:
            # File might have been deleted or is inaccessible
            return False
        
        # Check if we've seen this file before
        if filepath in self.file_sizes:
            previous_size, last_check_time = self.file_sizes[filepath]
            
            if current_size == previous_size:
                # Size hasn't changed - check if enough time has passed
                if time.time() - last_check_time >= self.stabilize_wait:
                    # File is stable
                    del self.file_sizes[filepath]  # Clean up tracking
                    return True
            else:
                # Size changed, update tracking
                self.file_sizes[filepath] = (current_size, time.time())
        else:
            # First time seeing this file
            self.file_sizes[filepath] = (current_size, time.time())
        
        return False
    
    def _process_dicom_file(self, dicom_path):
        """Convert DICOM file to NIfTI using dcm2niix."""
        try:
            print "Processing DICOM: %s" % dicom_path
            
            # Create temporary directory for dcm2niix output
            temp_dir = tempfile.mkdtemp(prefix='dcm2niix_')
            
            try:
                # Run dcm2niix
                cmd = [
                    'dcm2niix',
                    '-o', temp_dir,  # output directory
                    '-f', 'converted',  # output filename prefix
                    '-z', 'y',  # compress output
                    dicom_path
                ]
                
                print "Running: %s" % ' '.join(cmd)
                
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                stdout, stderr = proc.communicate()
                
                if proc.returncode == 0:
                    # Find the generated NIfTI file
                    nifti_files = glob.glob(os.path.join(temp_dir, '*.nii.gz'))
                    
                    if len(nifti_files) == 1:
                        # Move the file to our save location with sequential naming
                        temp_nifti = nifti_files[0]
                        
                        # Generate output filename
                        output_filename = os.path.join(
                            self.save_location,
                            'img-%05d.nii.gz' % self.save_volume_index
                        )
                        
                        shutil.move(temp_nifti, output_filename)
                        self.save_volume_index += 1
                        
                        print "Converted DICOM to: %s" % output_filename
                        
                        # Add to filename stack for processing
                        self.mutex.acquire()
                        self.filename_stack.append(output_filename)
                        self.mutex.release()
                        
                    elif len(nifti_files) == 0:
                        print "ERROR: dcm2niix did not produce any NIfTI files for %s" % dicom_path
                        print "stdout: %s" % stdout
                        print "stderr: %s" % stderr
                    else:
                        print "ERROR: dcm2niix produced multiple files for %s: %s" % (dicom_path, nifti_files)
                else:
                    print "ERROR: dcm2niix failed for %s (return code %d)" % (dicom_path, proc.returncode)
                    print "stdout: %s" % stdout
                    print "stderr: %s" % stderr
                
            finally:
                # Clean up temporary directory
                shutil.rmtree(temp_dir, ignore_errors=True)
                
        except Exception as e:
            print "ERROR processing DICOM file %s: %s" % (dicom_path, str(e))
        
        finally:
            # Remove from processing set
            self.processing_files.discard(dicom_path)
