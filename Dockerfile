# Use a base image with GCC, CMake, and Python
FROM gcc:latest

# Install CMake and Python
RUN apt-get update && apt-get install -y \
    cmake \
    python3 \
    python3-pip \
    python3-venv

# Set up Python Virtual Environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install Python scientific libraries
RUN pip install numpy matplotlib pandas

WORKDIR /app
COPY . .

# Build C++ Project
RUN mkdir build && cd build && cmake .. && make

# Default command: Run Visualization Pipeline
CMD ["/bin/bash", "-c", "./build/IsingModelSim --viz && python3 analysis/visualize.py"]