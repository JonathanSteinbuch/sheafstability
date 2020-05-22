
CPP_SRCS += \
./src/helpers/Helpers.cpp 

OBJS += \
./make/src/helpers/Helpers.o 

CPP_DEPS += \
./make/src/helpers/Helpers.d 


make/src/helpers/%.o: ./src/helpers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


