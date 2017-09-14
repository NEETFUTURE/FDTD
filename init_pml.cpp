#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#include <string>
#include "share.h"
#include "init_pml.h"
#include <math.h>


void InitPml
(
    //PML媒質内の係数
    double **CEzx,
    double **CEzx_x,
    double **CHyx,
    double **CHyx_x,
    double dz,
    double dt
)
{
    int M = 4;//導電率の分布を与える次数
    double Rn = 1.e-32f;
    double delta_max = -1.f * (double)(M + 1) * EPS0 * C * (double)logf(Rn) / 2.f / (double)PML / dz;
    double delta_E_X;
    double delta_H_X;
    *CEzx = (double* )malloc(sizeof(double)*PML);
    *CEzx_x = (double* )malloc(sizeof(double)*PML);
    *CHyx = (double* )malloc(sizeof(double)*PML);
    *CHyx_x = (double* )malloc(sizeof(double)*PML);

    for (int i=0; i < PML; i++)
    {
        delta_E_X = delta_max * (double)powf((double)(PML - i) / (double)PML, M);
        *(*CEzx+i)  = (1.f - delta_E_X * dt / 2.f / EPS0) / (1.f + delta_E_X * dt / 2.f / EPS0);
        *(*CEzx_x+i) = dt / EPS0 / (1.f + delta_E_X * dt / 2.f / EPS0) / dz;

        delta_H_X = delta_max * (double)powf(((double)(PML - i) - 1.f / 2.f) / (double)PML, M);
        *(*CHyx+i)   = (1.f - delta_H_X * dt / 2.f / EPS0) / (1.f + delta_H_X * dt / 2.f / EPS0);
        *(*CHyx_x+i) = dt / UMU / (1.f + delta_H_X * dt / 2.f / EPS0) / dz;
    }
}
