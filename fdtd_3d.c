#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 71
#define NY 71
#define NZ 71
#define NSTEP 600
#define UMU 1.257e-6
#define EPS0 8.854e-12
#define C 2.998e8
#define SIGMA 0.0
#define PI 3.14159265359
#define freq 0.6e9
#define STEP 1

#define TRGT_X 35
#define TRGT_Y 35
#define TRGT_Z 35

int main(int argc, char **argv)
{
    int n, i, j, k;

    double *hx = malloc(sizeof(double) * NX * NY * NZ);
    double *hy = malloc(sizeof(double) * NX * NY * NZ);
    double *hz = malloc(sizeof(double) * NX * NY * NZ);

    double *ex = malloc(sizeof(double) * NX * NY * NZ);
    double *ey = malloc(sizeof(double) * NX * NY * NZ);
    double *ez = malloc(sizeof(double) * NX * NY * NZ);

    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;

    char filename[20];
    char filename2[20];

    FILE *sample_fp;
    FILE *fp;

    t = 0.0;
    dt = dz / C/ 3.0;
    ec1 = (1.0 - SIGMA * dt / (2.0 * EPS0)) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    ec2 = dt / (EPS0 * dz) / (1.0 + SIGMA * dt / (2.0 * EPS0));
    hc = -dt / (dz * UMU);

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY; j++)
        {
            for (k = 0; k < NZ; k++)
            {
                ex[k*NX*NY + NX*j + i] = 0.0;
                ey[k*NX*NY + NX*j + i] = 0.0;
                ez[k*NX*NY + NX*j + i] = 0.0;
                hx[k*NX*NY + NX*j + i] = 0.0;
                hy[k*NX*NY + NX*j + i] = 0.0;
                hz[k*NX*NY + NX*j + i] = 0.0;
            }
        }
    }


    for (n = 0; n < NSTEP; n++)
    {

        // sprintf(filename, "data_3d/data_%04d.vtk", n);
        // fp = fopen(filename, "w");
        sprintf(filename2, "data_3d_raw/data_%04d.raw", n);
        sample_fp = fopen(filename2, "wb");

        // fprintf(fp, "# vtk DataFile Version 3.0\n");
        // fprintf(fp, "Example data 1\n");
        // fprintf(fp, "ASCII\n");
        // fprintf(fp, "DATASET STRUCTURED_POINTS\n");
        // fprintf(fp, "DIMENSIONS %d %d %d\n", NX+1, NY+1, NZ+1);
        // fprintf(fp, "SPACING 1 1 1\n");
        // fprintf(fp, "ORIGIN 0 0 0\n");
        // fprintf(fp, "CELL_DATA %d\n", NX*NY*NZ);
        // fprintf(fp, "SCALARS ez double\n");
        // fprintf(fp, "LOOKUP_TABLE default\n");

        fwrite(ez+TRGT_Z*NX*NY, sizeof(double), NX*NY, sample_fp);


        // ******************* 電界の計算 *********************
        if (t < 0.5 / freq)
        {
            ez[TRGT_Z*NX*NY + TRGT_Y*NX + TRGT_X] = ez[TRGT_Z*NX*NY + TRGT_Y*NX + TRGT_X] + pow(1.5*sin(2.0 * PI * freq * t), 4);
        }

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                for (k = 0; k < NZ; k++)
                {
                    if (j > 0 && k > 0)
                    {
                        ex[k*NX*NY + NX*j + i] = ec1 * ex[k*NX*NY + NX*j + i]
                                               + ec2 * (- hy[k*    NX*NY + NX*j     + i]
                                                        + hy[(k-1)*NX*NY + NX*j     + i]
                                                        + hz[k*    NX*NY + NX*j     + i]
                                                        - hz[k*    NX*NY + NX*(j-1) + i]);
                    }
                    if (i > 0 && k > 0)
                    {
                        ey[k*NX*NY + NX*j + i] = ec1 * ey[k*NX*NY + NX*j + i]
                                               + ec2 * (+ hx[k*    NX*NY + NX*j + i]
                                                        - hx[(k-1)*NX*NY + NX*j + i]
                                                        - hz[k*    NX*NY + NX*j + i]
                                                        + hz[k*    NX*NY + NX*j + i - 1]);
                    }
                    if (i > 0 && j > 0){
                        ez[k*NX*NY + NX*j + i] = ec1 * ez[k*NX*NY + NX*j + i]
                                               + ec2 * (- hx[k*NX*NY + NX*j     + i]
                                                        + hx[k*NX*NY + NX*(j-1) + i]
                                                        + hy[k*NX*NY + NX*j     + i]
                                                        - hy[k*NX*NY + NX*j     + i - 1]);
                    }
                    //fprintf(fp, "%9.7f ", ez[k*NX*NY + NX*j + i]);
                }
            }
            //fprintf(fp, "\n");
        }


        t = (t + dt / 2.0);

        // ******************* 磁界の計算 *********************
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                for (k = 0; k < NZ; k++)
                {
                    if (j < NY-1 && k < NZ-1){
                        hx[k*NX*NY + NX*j + i] = hx[k*NX*NY + NX*j + i]
                                               + hc * (+ ey[k*    NX*NY + NX*j     + i]
                                                       - ey[(k+1)*NX*NY + NX*j     + i]
                                                       - ez[k*    NX*NY + NX*j     + i]
                                                       + ez[k*    NX*NY + NX*(j+1) + i]);
                    }

                    if (i < NX-1 && k < NZ-1){
                        hy[k*NX*NY + NX*j + i] = hy[k*NX*NY + NX*j + i]
                                               + hc * (- ex[k*    NX*NY + NX*j + i]
                                                       + ex[(k+1)*NX*NY + NX*j + i]
                                                       + ez[k*    NX*NY + NX*j + i]
                                                       - ez[k*    NX*NY + NX*j + i + 1]);
                    }

                    if (j < NY-1 && i < NX-1){
                        hz[k*NX*NY + NX*j + i] = hz[k*NX*NY + NX*j + i]
                                               + hc * (+ ex[k*NX*NY + NX*j     + i]
                                                       - ex[k*NX*NY + NX*(j+1) + i]
                                                       - ey[k*NX*NY + NX*j     + i]
                                                       + ey[k*NX*NY + NX*j     + i + 1]);
                    }
                }
            }
        }

        t = (t + dt / 2.0);
        // fclose(fp);
        fclose(sample_fp);
    }


}
