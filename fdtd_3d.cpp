#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>

#define NX 71
#define NY 71
#define NZ 71
#define NSTEP 200
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

    double *hx = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *hy = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *hz = (double*)malloc(sizeof(double) * NX * NY * NZ);

    double *ex = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *ey = (double*)malloc(sizeof(double) * NX * NY * NZ);
    double *ez = (double*)malloc(sizeof(double) * NX * NY * NZ);

    double dt, ec1, ec2, hc;
    double t;
    double dz = 1.0e-2;

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


        // ******************* 電界の計算 *********************
        
        if (t < 4.0 / freq)
        for(i=0;i<NZ;i++){
            {
                ez[i*NX*NY + TRGT_Y*NX + TRGT_X] = ez[i*NX*NY + TRGT_Y*NX + TRGT_X] + pow(sin(2.0 * PI * freq * t), 4);
            }
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
                    // fprintf(fp, "%9.7f ", ez[k*NX*NY + NX*j + i]);
                    double* pixel = static_cast<double*>(imageData->GetScalarPointer(i,j,k));
                    pixel[0] = ez[k*NX*NY + NX*j + i];
                }
            }
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
