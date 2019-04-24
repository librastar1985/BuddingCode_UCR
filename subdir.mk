################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11

current_dir := $(shell pwd)
LIBS:=  -lpugixml -L/$(current_dir)/pugixml/lib64
#-lgsl -lgslcblas

ILIBS_cuda8 = -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/8.0/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AreaTriangles.cu \
../BendingTriangles.cu \
../MemRepulsionSprings.cu\
../LinearSprings.cu \
../LJSprings.cu \
../NodeAdvance.cu \
../AreaTrianglesEnergy.cu \
../BendingTrianglesEnergy.cu \
../LinearSpringsEnergy.cu \
../LJEnergy.cu \
../MemRepulsionEnergy.cu \
../System.cu \
../Edgeswap_test.cpp \
../SystemBuilder.cpp \
../Storage.cu \
../main.cpp


# this is a variable
OBJS += \
./AreaTriangles.o \
./BendingTriangles.o \
./MemRepulsionSprings.o\
./LinearSprings.o \
./LJSprings.o \
./NodeAdvance.o \
./AreaTrianglesEnergy.o \
./BendingTrianglesEnergy.o \
./LinearSpringsEnergy.o \
./LJEnergy.o \
./MemRepulsionEnergy.o \
./System.o \
./Edgeswap_test.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./AreaTriangles.d \
./BendingTriangles.d \
./MemRepulsionSprings.d\
./LinearSprings.d \
./LJSprings.d \
./NodeAdvance.d \
./AreaTrianglesEnergy.d \
./BendingTrianglesEnergy.d \
./LinearSpringsEnergy.d \
./LJEnergy.d \
./MemRepulsionEnergy.d \
./System.d \
./Edgeswap_test.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBS_cuda8) $(LIBS) -o $@ -c $^ 

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS_cuda8) -dc -o $@ $^
