# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./src/DensePolyMatrix.cpp\
./src/Helpers.cpp \
./src/Poly.cpp \
./src/PolyRing.cpp \
./src/SparsePolyMatrix.cpp \
./src/SparseScalarMatrix.cpp \
./src/stability.cpp 

OBJS += \
./make/src/DensePolyMatrix.o \
./make/src/Helpers.o \
./make/src/Poly.o \
./make/src/PolyRing.o \
./make/src/SparsePolyMatrix.o \
./make/src/SparseScalarMatrix.o \
./make/src/stability.o 

CPP_DEPS += \
./make/src/DensePolyMatrix.d \
./make/src/Helpers.d \
./make/src/Poly.d \
./make/src/PolyRing.d \
./make/src/SparsePolyMatrix.d \
./make/src/SparseScalarMatrix.d \
./make/src/stability.d 


# Each subdirectory must supply rules for building sources it contributes
make/src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


