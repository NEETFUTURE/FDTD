#include "share.h"
#include "init_pml.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>

const int NX = 150;
const int NY = 70;
const int NZ = 10;
const int PML = 16;
const int NSTEP = 400;
const double UMU = 1.257e-6;
const double EPS0 = 8.854e-12;
const double C = 2.998e8;
const double SIGMA = 0.0;
const double PI = 3.14159265359;
const double freq = 0.6e9;
const int STEP = 1;

const int TRGT_X = 50;
const int TRGT_Y = 35;
const int TRGT_Z = 35;


int main(int argc, char **argv)
{
    int n, i, j, k;

    double *hx = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *hy = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *hz = (double*)malloc(sizeof(double) * NX * NY * NZ);

    double *ex = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *ey = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *ez = (double*)malloc(sizeof(double) * NX * NY * NZ);

    //PML層内の電界・磁界
    double r_Ezx;
    double *Ezx = (double*)malloc(sizeof(double) * PML * NX * NY);
    double r_Hyx;
    double *Hyx = (double*)malloc(sizeof(double) * PML * NX * NY);

    //PML媒質内の係数
    double *CEzx;
    double *CEzx_x;

    double *CHyx;
    double *CHyx_x;

    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;
	double omega = 2*PI*freq;

    char filename[40];

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
    memset(Ezx, 0, sizeof(double) * PML * NX * NY);
    memset(Hyx, 0, sizeof(double) * PML * NX * NY);
    InitPml
    (
        //PML媒質内の係数
        &CEzx,
        &CEzx_x,
        &CHyx,
        &CHyx_x,
        dz,
        dt
    );
    printf("Init_pml\n");

    vtkSmartPointer<vtkImageData> imageData =
    vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(NX+1, NY+1, NZ+1);
#if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToDouble();
#else
    imageData->AllocateScalars(VTK_DOUBLE, 1);
#endif
    int* dims = imageData->GetDimensions();

    for (n = 0; n < NSTEP; n++)
    {

        printf("Done n = %5d/%5d\n", n, NSTEP);
        sprintf(filename, "../data/data_%04d.vti", n);


		// ******************* 波源 ******************************
		if (t < 0.5 / freq)
		for (k = 0; k<NZ - 1; k++) {
			for (j = 0; j<NY; j++) {
				ez[k*NX*NY + j*NX + TRGT_X] = ez[k*NX*NY + j*NX + TRGT_X] + powf(sin(omega * t), 4.0);
			}
		}

		// ******************* 完全導体 **************************
		//int Xr = (int)(NX*2/3);
		//int Yr = (int)(NY/2);
		//int Rr = 5;

		//for (i = Xr - Rr; i <= Xr + Rr; i++)
		//{
		//	for (j = Yr - Rr; j <= Yr + Rr; j++)
		//	{
		//		if (pow(i - Xr, 2.0) + pow(j - Yr, 2.0) <= Rr*Rr)
		//		{
		//			for (k = 0; k < NZ - 1; k++)
		//			{
		//				ey[k*NX*NY + NX*j + i] = 0.0;
		//				ex[k*NX*NY + NX*j + i] = 0.0;
		//				ez[k*NX*NY + NX*j + i] = 0.0;
		//			}
		//		}
		//	}
		//}

		// ******************* 電界PMLの計算 *********************

		for (i = 1; i < PML; i++)
		{
			for (j = 1; j < NY; j++)
			{
				for (k = 0; k < NZ - 1; k++)
				{
					ez[k*NX*NY + NX*j + i] = CEzx[i] * ez[k*NX*NY + NX*j + i]
					                       + CEzx_x[i] * (hy[k*NX*NY + NX*j + i] - hy[k*NX*NY + NX*j + i - 1]);
					double* pixel = static_cast<double*>(imageData->GetScalarPointer(i, j, k));
					pixel[0] = ez[k*NX*NY + NX*j + i];

				}
			}
		}
        // ******************* 電界の計算 *********************


        for (i = PML; i < NX-1; i++)
        {
            for (j = 1; j < NY; j++)
            {
                for (k = 0; k < NZ-1; k++)
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

                    double* pixel = static_cast<double*>(imageData->GetScalarPointer(i,j,k));
                    pixel[0] = ez[k*NX*NY + NX*j + i];
                }
            }
        }


        t = (t + dt / 2.0);

        // ******************* 磁界の計算 *********************
        for (i = PML; i < NX; i++)
        {
            for (j = 1; j < NY; j++)
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
                                               + hc * (+ ex[k*    NX*NY + NX*j + i]
												       - ex[(k+1)*NX*NY + NX*j + i]
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
        // printf("磁界計算\n");
        // ******************* 磁界PMLの計算 *********************

        for (i = 1; i < PML; i++)
        {
            for (j = 1; j < NY; j++)
            {
                for (k = 0; k < NZ; k++)
                {
                    if (i < NX-1 && k < NZ-1){
                        hy[k*NX*NY + NX*j + i] = CHyx[i] * hy[k*NX*NY + NX*j + i]
                                               + CHyx_x[i] * (-ez[k*NX*NY + NX*j + i + 1]
                                                              + ez[k*NX*NY + NX*j + i]);
                    }
                }
            }
        }

        // printf("磁界PML計算\n");
        t = (t + dt / 2.0);



        vtkSmartPointer<vtkXMLImageDataWriter> writer =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();
        writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
#else
        writer->SetInputData(imageData);
#endif
        writer->Write();
    }
}
