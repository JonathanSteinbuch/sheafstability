################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/helpers/Helpers.cpp 

OBJS += \
./src/helpers/Helpers.o 

CPP_DEPS += \
./src/helpers/Helpers.d 


# Each subdirectory must supply rules for building sources it contributes
src/helpers/%.o: ../src/helpers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -I/home/math/jsteinbu/lela/include -O0 -fopenmp -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


