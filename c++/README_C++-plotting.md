# ...existing code...

## Installing matplotlib-cpp

1. Install Python and matplotlib:
   ```sh
   sudo apt-get update
   sudo apt-get install python3 python3-pip
   pip3 install matplotlib
   ```

2. Clone the `matplotlib-cpp` repository:
   ```sh
   git clone https://github.com/lava/matplotlib-cpp.git
   ```

3. Install Python development headers:
   ```sh
   sudo apt-get install python3-dev
   ```

4. Update the `CMakeLists.txt` file to include the `matplotlib-cpp` directory and link against Python libraries.

5. Build the project:
   ```sh
   mkdir build
   cd build
   cmake ..
   make
   ```

# ...existing code...
