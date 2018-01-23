# Matrix

### Description

This is a basic implementation of a ground-up Matrix class, with operations supporting accessing elements at indices (i,j), transpose, and matrix multiplication. The class is templated to accept any numeric data type, and holds he data in an underlying 1D STL vector.

### Testing

Included with the matrix.h file, one may also used the provided CMakeList.txt file to run a series of tests using the Googletest framework in an exectuable file called **testAll**. These tests include basic instantiation, transpose testing, matrix multiplication, and performance benchmarking against vanilla, naive versions of transpose and matrix multiplication
