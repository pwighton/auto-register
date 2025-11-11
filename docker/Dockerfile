FROM pwighton/freesurfer:8.1.0

WORKDIR /opt/auto-register

# Copy FreeSurfer license file
COPY docker/.license /usr/local/freesurfer/8.1.0-1/.license

# Install Python 2.7, pip, unzip, and wget
# Using python2 specifically to avoid disturbing any existing python environment
RUN dnf -y install python2 python2-pip unzip wget && \
    dnf clean all

# Download and install dcm2niix pre-built binary
RUN wget https://github.com/rordenlab/dcm2niix/releases/download/v1.0.20250506/dcm2niix_lnx.zip -O /tmp/dcm2niix.zip && \
    unzip /tmp/dcm2niix.zip -d /tmp/dcm2niix && \
    mv /tmp/dcm2niix/dcm2niix /usr/local/bin/dcm2niix && \
    chmod +x /usr/local/bin/dcm2niix && \
    rm -rf /tmp/dcm2niix.zip /tmp/dcm2niix

# Install Python 2 dependencies
# Using python2 -m pip to ensure we're using Python 2's pip
COPY requirements.txt /opt/auto-register/requirements.txt
RUN python2 -m pip install --no-cache-dir -r requirements.txt

# Copy auto-register source code
COPY src/ /opt/auto-register/src/

# Add auto-register to PATH
ENV PATH="/opt/auto-register/src:${PATH}"

# Set default command
CMD ["/bin/bash"]
