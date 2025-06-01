To build the program you need a C++ Compiler with OpenMP support and CMake and a Fortran compiler. We have test this on Unix based systems (including WSL on Windows). To compile the program 
    mkdir build
    cd build
    cmake ..
    make

This generates an executable vp-solver that can then be run as follows
    OMP_NUM_THREADS=4 ./vp-solver

You can plot the time evolution of the electric field easily using gnuplot
    gnuplot
    gnuplot> p "evolution.data" u 1:2

    

