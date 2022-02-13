#include <stdlib.h>
#include "simulator.hpp"
#include <iostream>

int main(int argc, char** argv) {
    // Size from command line if specified
    unsigned grid = 256;
    if (argc > 1) {
        auto tmp = atoi(argv[1]);
        if(tmp <= 0) {
            std::cerr << "Grid should be larger than 0\n";
            return EXIT_FAILURE;
        }
        grid = static_cast<decltype(grid)>(tmp);
    }
    std::cout << "The size of the grid is: " << grid << "\n";
    Simulator s{grid};
    s.run(4.5, 1.0, 100.0);


    // double* u = new double[grid * (grid + 1)];
    // double* un = new double[grid * (grid + 1)];
    // double* uc = new double[grid * grid];
    // double* v = new double[(grid + 1) * grid];
    // double* vn = new double[(grid + 1) * grid];
    // double* vc = new double[grid * grid];
    // double* p = new double[(grid + 1) * (grid + 1)];
    // double* pn = new double[(grid + 1) * (grid + 1)];
    // double* pc = new double[grid * grid];
    // double* m = new double[(grid + 1) * (grid + 1)];

    // int i, j, step;
    // double dx, dy, dt, tau, delta, error, Re;
    // step = 1;
    // dx = 1.0 / (grid - 1);
    // dy = 1.0 / (grid - 1);
    // dt = 0.001 / ((grid / 128 * 2) * (grid / 128 * 2));
    // delta = 4.5;
    // error = 1.0;
    // Re = 100.0;

    // // Initializing u
    // for (i = 0; i <= (grid - 1); i++) {
    //     for (j = 0; j <= (grid); j++) {
    //         u[(i) * (grid + 1) + j] = 0.0;
    //         u[(i) * (grid + 1) + grid] = 1.0;
    //         u[(i) * (grid + 1) + grid - 1] = 1.0;
    //     }
    // }

    // // Initializing v
    // for (i = 0; i <= (grid); i++) {
    //     for (j = 0; j <= (grid - 1); j++) { v[(i)*grid + j] = 0.0; }
    // }

    // // Initializing p
    // for (i = 0; i <= (grid); i++) {
    //     for (j = 0; j <= (grid); j++) { p[(i) * (grid + 1) + j] = 1.0; }
    // }

    // while (error > 0.00000001) {
    //     // Solve u-momentum equation
    //     for (i = 1; i <= (grid - 2); i++) {
    //         for (j = 1; j <= (grid - 1); j++) {
    //             un[(i) * (grid + 1) + j] =
    //                 u[(i) * (grid + 1) + j]
    //                 - dt
    //                       * ((u[(i + 1) * (grid + 1) + j] * u[(i + 1) * (grid + 1) + j]
    //                              - u[(i - 1) * (grid + 1) + j] * u[(i - 1) * (grid + 1) + j])
    //                               / 2.0 / dx
    //                           + 0.25
    //                                 * ((u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j + 1])
    //                                         * (v[(i)*grid + j] + v[(i + 1) * grid + j])
    //                                     - (u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j - 1])
    //                                           * (v[(i + 1) * grid + j - 1] + v[(i)*grid + j - 1]))
    //                                 / dy)
    //                 - dt / dx * (p[(i + 1) * (grid + 1) + j] - p[(i) * (grid + 1) + j])
    //                 + dt * 1.0 / Re
    //                       * ((u[(i + 1) * (grid + 1) + j] - 2.0 * u[(i) * (grid + 1) + j]
    //                              + u[(i - 1) * (grid + 1) + j])
    //                               / dx / dx
    //                           + (u[(i) * (grid + 1) + j + 1] - 2.0 * u[(i) * (grid + 1) + j]
    //                                 + u[(i) * (grid + 1) + j - 1])
    //                                 / dy / dy);
    //         }
    //     }

    //     // Boundary conditions
    //     for (j = 1; j <= (grid - 1); j++) {
    //         un[(0) * (grid + 1) + j] = 0.0;
    //         un[(grid - 1) * (grid + 1) + j] = 0.0;
    //     }

    //     for (i = 0; i <= (grid - 1); i++) {
    //         un[(i) * (grid + 1) + 0] = -un[(i) * (grid + 1) + 1];
    //         un[(i) * (grid + 1) + grid] = 2 - un[(i) * (grid + 1) + grid - 1];
    //     }

    //     // Solves v-momentum
    //     for (i = 1; i <= (grid - 1); i++) {
    //         for (j = 1; j <= (grid - 2); j++) {
    //             vn[(i)*grid + j] =
    //                 v[(i)*grid + j]
    //                 - dt
    //                       * (0.25
    //                               * ((u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j + 1])
    //                                       * (v[(i)*grid + j] + v[(i + 1) * grid + j])
    //                                   - (u[(i - 1) * (grid + 1) + j]
    //                                         + u[(i - 1) * (grid + 1) + j + 1])
    //                                         * (v[(i)*grid + j] + v[(i - 1) * grid + j]))
    //                               / dx
    //                           + (v[(i)*grid + j + 1] * v[(i)*grid + j + 1]
    //                                 - v[(i)*grid + j - 1] * v[(i)*grid + j - 1])
    //                                 / 2.0 / dy)
    //                 - dt / dy * (p[(i) * (grid + 1) + j + 1] - p[(i) * (grid + 1) + j])
    //                 + dt * 1.0 / Re
    //                       * ((v[(i + 1) * grid + j] - 2.0 * v[(i)*grid + j] + v[(i - 1) * grid + j])
    //                               / dx / dx
    //                           + (v[(i)*grid + j + 1] - 2.0 * v[(i)*grid + j] + v[(i)*grid + j - 1])
    //                                 / dy / dy);
    //         }
    //     }

    //     // Boundary conditions
    //     for (j = 1; j <= (grid - 2); j++) {
    //         vn[(0) * grid + j] = -vn[(1) * grid + j];
    //         vn[(grid)*grid + j] = -vn[(grid - 1) * grid + j];
    //     }

    //     for (i = 0; i <= (grid); i++) {
    //         vn[(i)*grid + 0] = 0.0;
    //         vn[(i)*grid + grid - 1] = 0.0;
    //     }

    //     // Solves continuity equation
    //     for (i = 1; i <= (grid - 1); i++) {
    //         for (j = 1; j <= (grid - 1); j++) {
    //             pn[(i) * (grid + 1) + j] =
    //                 p[(i) * (grid + 1) + j]
    //                 - dt * delta
    //                       * ((un[(i) * (grid + 1) + j] - un[(i - 1) * (grid + 1) + j]) / dx
    //                           + (vn[(i)*grid + j] - vn[(i)*grid + j - 1]) / dy);
    //         }
    //     }

    //     // Boundary conditions
    //     for (i = 1; i <= (grid - 1); i++) {
    //         pn[(i) * (grid + 1) + 0] = pn[(i) * (grid + 1) + 1];
    //         pn[(i) * (grid + 1) + grid] = pn[(i) * (grid + 1) + grid - 1];
    //     }

    //     for (j = 0; j <= (grid); j++) {
    //         pn[(0) * (grid + 1) + j] = pn[(1) * (grid + 1) + j];
    //         pn[(grid) * (grid + 1) + j] = pn[(grid - 1) * (grid + 1) + j];
    //     }

    //     // Displaying error
    //     error = 0.0;

    //     for (i = 1; i <= (grid - 1); i++) {
    //         for (j = 1; j <= (grid - 1); j++) {
    //             m[(i) * (grid + 1) + j] =
    //                 ((un[(i) * (grid + 1) + j] - un[(i - 1) * (grid + 1) + j]) / dx
    //                     + (vn[(i)*grid + j] - vn[(i)*grid + j - 1]) / dy);
    //             error = error + fabs(m[(i) * (grid + 1) + j]);
    //         }
    //     }

    //     // if (step%1000 ==1)
    //     { printf("Error is %5.8lf for the step %d\n", error, step); }

    //     // Iterating u
    //     for (i = 0; i <= (grid - 1); i++) {
    //         for (j = 0; j <= (grid); j++) { u[(i) * (grid + 1) + j] = un[(i) * (grid + 1) + j]; }
    //     }

    //     // Iterating v
    //     for (i = 0; i <= (grid); i++) {
    //         for (j = 0; j <= (grid - 1); j++) { v[(i)*grid + j] = vn[(i)*grid + j]; }
    //     }

    //     // Iterating p
    //     for (i = 0; i <= (grid); i++) {
    //         for (j = 0; j <= (grid); j++) { p[(i) * (grid + 1) + j] = pn[(i) * (grid + 1) + j]; }
    //     }

    //     step = step + 1;
    // }

    // for (i = 0; i <= (grid - 1); i++) {
    //     for (j = 0; j <= (grid - 1); j++) {
    //         uc[(i)*grid + j] = 0.5 * (u[(i) * (grid + 1) + j] + u[(i) * (grid + 1) + j + 1]);
    //         vc[(i)*grid + j] = 0.5 * (v[(i)*grid + j] + v[(i + 1) * grid + j]);
    //         pc[(i)*grid + j] = 0.25
    //                            * (p[(i) * (grid + 1) + j] + p[(i + 1) * (grid + 1) + j]
    //                                + p[(i) * (grid + 1) + j + 1] + p[(i + 1) * (grid + 1) + j + 1]);
    //     }
    // }

    // // OUTPUT DATA
    // FILE *fout2, *fout3;
    // fout2 = fopen("UVP.plt", "w+t");
    // fout3 = fopen("Central_U.plt", "w+t");

    // if (fout2 == NULL) {
    //     printf("\nERROR when opening file\n");
    //     fclose(fout2);
    // }

    // else {
    //     fprintf(fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
    //     fprintf(fout2, "ZONE  F=POINT\n");
    //     fprintf(fout2, "I=%d, J=%d\n", grid, grid);

    //     for (j = 0; j < (grid); j++) {
    //         for (i = 0; i < (grid); i++) {
    //             double xpos, ypos;
    //             xpos = i * dx;
    //             ypos = j * dy;

    //             fprintf(fout2,
    //                 "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n",
    //                 xpos,
    //                 ypos,
    //                 uc[(i)*grid + j],
    //                 vc[(i)*grid + j],
    //                 pc[(i)*grid + j]);
    //         }
    //     }
    // }

    // fclose(fout2);

    // // CENTRAL --U
    // fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
    // fprintf(fout3, "ZONE F=POINT\n");
    // fprintf(fout3, "I=%d\n", grid);

    // for (j = 0; j < grid; j++) {
    //     double ypos;
    //     ypos = (double)j * dy;

    //     fprintf(fout3,
    //         "%5.8lf\t%5.8lf\n",
    //         (uc[(grid / 2) * grid + j] + uc[((grid / 2) + 1) * grid + j]) / (2.),
    //         ypos);
    // }

    // delete[] u;
    // delete[] un;
    // delete[] uc;
    // delete[] v;
    // delete[] vn;
    // delete[] vc;
    // delete[] p;
    // delete[] pn;
    // delete[] pc;
    // delete[] m;
}
