#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define N 2048
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int cannon(int argc, char* argv[]) {
    int dim, numtasks, rank, source, dest, outbuf, i, tag = 1, inbuf[4] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL }, nbrs[4],
                                                 dims[2], periods[2] = { 0, 0 }, reorder = 0, coords[2];
    MPI_Request reqs[8];
    MPI_Status stats[8];
    MPI_Comm cartcomm;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout<<"The rank is: "<<rank<<std::endl;
        long root = sqrt(numtasks);
        if (root * root == numtasks) {
            dim = root;
            printf("Number of tasks is a square number: %i \n", numtasks);
        } else {
            std::cout<<rank<<std::endl;
            printf("Number of tasks isn't a square number: %i \n", numtasks);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Use MPI_Dims_create to figure out how many processes in each direction
    dims[0] = dims[1] = dim;
    MPI_Dims_create(numtasks, 2, dims);
    // Create cartesian communicator for 2D, dims[0]*dims[1] processes,
    // without periodicity and reordering
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm);
    // Get my rank in the new communicator
    MPI_Comm_rank(cartcomm, &rank);
    // Get my coordinates coords[0] and coords[1]
    MPI_Cart_coords(cartcomm, rank, 2, coords);
    // Get my neighbours in dimension 0
    MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
    // Get my neighbours in dirmension 1
    MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);
    //printf("rank= %d coords= %d %d  neighbors(u,d,l,r)= %d %d %d %d\n", rank, coords[0], coords[1], nbrs[UP], nbrs[DOWN], nbrs[LEFT], nbrs[RIGHT]);
    outbuf = rank;

    
    //MPI_Waitall(8, reqs, stats);
    MPI_Finalize();
    return 0;
}
