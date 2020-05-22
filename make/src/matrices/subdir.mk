################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/matrices/DensePolyMatrix.cpp \
../src/matrices/SparsePolyMatrix.cpp \
../src/matrices/SparseScalarMatrix.cpp 

OBJS += \
./src/matrices/DensePolyMatrix.o \
./src/matrices/SparsePolyMatrix.o \
./src/matrices/SparseScalarMatrix.o 

CPP_DEPS += \
./src/matrices/DensePolyMatrix.d \
./src/matrices/SparsePolyMatrix.d \
./src/matrices/SparseScalarMatrix.d 


# Each subdirectory must supply rules for building sources it contributes
src/matrices/%.o: ../src/matrices/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -I/home/math/jsteinbu/lela/include -O0 -fopenmp -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


