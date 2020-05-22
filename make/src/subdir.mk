
CPP_SRCS += \
./src/stability.cpp 

OBJS += \
./make/src/stability.o 

CPP_DEPS += \
./make/src/stability.d 


make/src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


