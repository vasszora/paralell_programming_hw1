#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define N 200
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int cannon(int argc, char* argv[]) {
    typedef float real;
    int dim, numtasks, rank, dims[2], periods[2] = { 1, 1 }, reorder = 0, coords[2], size_per_process, nbrs[4];
    std::vector<real> a(N*N);
    std::vector<real> b(N*N);
    std::vector<real> c(N*N);
    std::vector<real> buffer, temp;
    MPI_Comm cartcomm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        long root = sqrt(numtasks);
        if (root * root == numtasks) {
            printf("Number of tasks is a square number: %i \n", numtasks);
            dim = root;
            size_per_process = N / dim;
        } else {
            printf("Number of tasks isn't a square number: %i \n", numtasks);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size_per_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Use MPI_Dims_create to figure out how many processes in each direction
    dims[0] = dims[1] = dim;
    MPI_Dims_create(numtasks, 2, dims);
    // Create cartesian communicator for 2D, dims[0]*dims[1] processes,
    // without and reordering and with periodicity
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
    // Get my rank in the new communicator
    MPI_Comm_rank(cartcomm, &rank);
    // Get my coordinates coords[0] and coords[1]
    MPI_Cart_coords(cartcomm, rank, 2, coords);
 
    //Fill matrices with random values
    for (int i = coords[0]*size_per_process; i < (coords[0] + 1)*size_per_process; i++) {
        for (int j = coords[1]*size_per_process; j < (coords[1] + 1)*size_per_process; j++) {
            a[i*size_per_process+j] = (real)std::rand()/(real)RAND_MAX;
            b[i*size_per_process+j] = (real)std::rand()/(real)RAND_MAX;
        }
    }

    //initial shifting
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[LEFT],  &nbrs[RIGHT]); 
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[UP], &nbrs[DOWN]); 
    
    buffer.resize(size_per_process*size_per_process);
    //Calculate multiplication for subarray
    for (int shift = 0; shift < dim; shift++){
        for (int i = coords[0]*size_per_process; i < (coords[0] + 1)*size_per_process; i++) {
            for (int j = coords[1]*size_per_process; j < (coords[1] + 1)*size_per_process; j++)
                for (int k = 0; k < size_per_process; k++) {
                    c[i*size_per_process+j] += b[i*size_per_process+k] * a[k*size_per_process+j];
            }
        }
        if ( shift == dims[0]-1 ) break;
        //Shift matrices
        MPI_Sendrecv(&a, size_per_process*size_per_process, MPI_DOUBLE, nbrs[LEFT], 1, &buffer, size_per_process*size_per_process, MPI_DOUBLE, nbrs[RIGHT], 1, cartcomm, &status); 
        temp = buffer; buffer = a; b= temp;
        MPI_Sendrecv(&b, size_per_process*size_per_process, MPI_DOUBLE, nbrs[UP], 1, &buffer, size_per_process*size_per_process, MPI_DOUBLE, 1, nbrs[DOWN], cartcomm, &status); 
        temp = buffer; buffer = b; b = temp;
    }
    
    MPI_Finalize();
    return 0;
}
