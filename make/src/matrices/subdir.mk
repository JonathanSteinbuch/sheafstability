
CPP_SRCS += \
./src/matrices/DensePolyMatrix.cpp \
./src/matrices/SparsePolyMatrix.cpp \
./src/matrices/SparseScalarMatrix.cpp 

OBJS += \
./make/src/matrices/DensePolyMatrix.o \
./make/src/matrices/SparsePolyMatrix.o \
./make/src/matrices/SparseScalarMatrix.o 

CPP_DEPS += \
./make/src/matrices/DensePolyMatrix.d \
./make/src/matrices/SparsePolyMatrix.d \
./make/src/matrices/SparseScalarMatrix.d 


make/src/matrices/%.o: ./src/matrices/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


