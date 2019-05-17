# the base ubuntu LTS image
FROM ubuntu:18.04

# Install Ubuntu package dependencies
RUN apt-get update && apt-get install -y \
    python3-numpy \
    python3-dev \
    python3-scipy \
    python3-nose \
    python3-pip \
    gfortran \
    make \
    gmsh \
    libatlas-base-dev \
    libblas-dev \
    liblapack-dev \
    libsuitesparse-dev \
    ssh

# Run the pip install for matplotlib
RUN pip3 install matplotlib

# Clean up the apt-getting
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Add the NumBAT source code -> Will be overwritten when used with mounted volumes! (which is good)
COPY ./ /home/NumBAT/

# Compile the Fortran code
WORKDIR /home/NumBAT/backend/fortran/
RUN make

# Add the backend files to the python path
ENV PYTHONPATH "${PYTHONPATH}:/home/NumBAT/backend/"

# Add the matplotlib stylefile to the correct location
WORKDIR /root/.config/matplotlib
RUN mkdir stylelib && cp /home/NumBAT/backend/NumBATstyle.mplstyle stylelib/

# Run the tests
WORKDIR /home/NumBAT/tests/
RUN nosetests3

# Change the working directory to final spot
WORKDIR /home/NumBAT/

# Finish with the bash
CMD ["/bin/bash"]