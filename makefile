TARGET1 = fdtd1d
TARGET2 = fdtd_tm

all: $(TARGET1) $(TARGET2)
$(TARGET1): fdtd1d.f03
	gfortran -o $(TARGET1) -g -std=f2003 -Wall -fbounds-check fdtd1d.f03
$(TARGET2): fdtd_tm.f03
	gfortran -o $(TARGET2) -g -std=f2003 -Wall -fbounds-check fdtd_tm.f03

