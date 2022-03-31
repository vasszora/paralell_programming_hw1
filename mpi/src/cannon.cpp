#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define N 4 
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int cannon(int argc, char* argv[]) {
    typedef float real;
    int dim, numtasks, rank, source, dest, outbuf, i, tag = 1, inbuf[4] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL }, nbrs[4],
                                                 dims[2], periods[2] = { 1, 1 }, reorder = 0, coords[2];
    int size_per_process;
    std::vector<real> a(N*N);
    std::vector<real> b(N*N);
    std::vector<real> c(N*N);
    MPI_Request reqs[8];
    MPI_Status stats[8];
    MPI_Comm cartcomm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

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
    // std::cout<<coords[0]<<std::endl;
    // std::cout<<coords[1]<<std::endl;
    // std::cout<<"------------------"<<std::endl;
    // Get my neighbours in dimension 0
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    // Get my neighbours in dirmension 1
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
    //printf("rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n", rank, coords[0], coords[1], nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);

    //Fill matrices with random values
    for (int i = coords[0]*size_per_process; i < (coords[0] + 1)*size_per_process; i++) {
        for (int j = coords[1]*size_per_process; j < (coords[1] + 1)*size_per_process; j++) {
            a[i*N+j] = (real)std::rand()/(real)RAND_MAX;
            b[i*N+j] = (real)std::rand()/(real)RAND_MAX;
        }
    }
    MPI_Barrier(cartcomm);

    //Calculate multiplication

    //MPI_Waitall(8, reqs, stats);
    MPI_Finalize();
    return 0;
}
