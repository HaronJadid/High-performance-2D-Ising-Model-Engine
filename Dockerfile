# Use a base image with GCC and CMake installed
FROM gcc:latest

# Install CMake
RUN apt-get update && apt-get install -y cmake

# Set the working directory
WORKDIR /app

# Copy the source code
COPY . .

# Create a build directory
RUN mkdir build

# Build the project
WORKDIR /app/build
RUN cmake .. && make

# Command to run the application (optional, adjust as needed)
CMD ["./IsingModelSim"]