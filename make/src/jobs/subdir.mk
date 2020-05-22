################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/jobs/JacobiTaylorJob.cpp \
../src/jobs/JobType.cpp \
../src/jobs/KernelJob.cpp \
../src/jobs/PowersJob.cpp \
../src/jobs/ReaderAndCaller.cpp \
../src/jobs/SemistabilityJob.cpp 

OBJS += \
./src/jobs/JacobiTaylorJob.o \
./src/jobs/JobType.o \
./src/jobs/KernelJob.o \
./src/jobs/PowersJob.o \
./src/jobs/ReaderAndCaller.o \
./src/jobs/SemistabilityJob.o 

CPP_DEPS += \
./src/jobs/JacobiTaylorJob.d \
./src/jobs/JobType.d \
./src/jobs/KernelJob.d \
./src/jobs/PowersJob.d \
./src/jobs/ReaderAndCaller.d \
./src/jobs/SemistabilityJob.d 


# Each subdirectory must supply rules for building sources it contributes
src/jobs/%.o: ../src/jobs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -I/home/math/jsteinbu/lela/include -O0 -fopenmp -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


