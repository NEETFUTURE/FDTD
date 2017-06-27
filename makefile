TARGET1 = fdtd_1d
TARGET2 = fdtd_tm
TARGET3 = fdtd_3d

all: $(TARGET1) $(TARGET2) $(TARGET3)
$(TARGET1): fdtd_1d.f03
	gfortran -o $(TARGET1) -g -std=f2003 -Wall -fbounds-check fdtd_1d.f03
$(TARGET2): fdtd_tm.f03
	gfortran -o $(TARGET2) -g -std=f2003 -Wall -fbounds-check fdtd_tm.f03
$(TARGET3): fdtd_3d.c
	gcc fdtd_3d.c -o $(TARGET3) -lm
