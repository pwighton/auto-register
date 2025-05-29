class ImageManager(object):
    """ImageManager that supports both TCP socket and directory monitoring modes."""
    
    def __init__(self, args):
        """Initialize receiver based on input mode."""
        self.input_mode = getattr(args, 'input_mode', 'socket')
        
        if self.input_mode == 'socket':
            # Import and initialize the original ImageReceiver
            from image_receiver import ImageReceiver
            self.receiver = ImageReceiver(args)
        elif self.input_mode == 'directory':
            # Use directory monitoring
            from directory_monitor import DirectoryMonitor
            if not hasattr(args, 'watch_directory') or args.watch_directory is None:
                raise ValueError("--watch-directory must be specified when using directory input mode")
            self.receiver = DirectoryMonitor(args)
        else:
            raise ValueError("Invalid input_mode: %s. Must be 'socket' or 'directory'" % self.input_mode)
    
    def start(self):
        """Start the receiver."""
        self.receiver.start()
    
    def stop(self):
        """Stop the receiver."""
        self.receiver.stop()
    
    def get_next_filename(self):
        """Get the next available filename."""
        return self.receiver.get_next_filename()
    
    def is_running(self):
        """Check if receiver is running."""
        return self.receiver.is_running()
