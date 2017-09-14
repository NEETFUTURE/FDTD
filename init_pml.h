#ifndef OUTPUT_H_
#define OUTPUT_H_

void InitPml
(
    //PML媒質内の係数
    double **CEzx,
    double **CEzx_x,
    double **CHyx,
    double **CHyx_x,
    double dz,
    double dt
);

#endif
