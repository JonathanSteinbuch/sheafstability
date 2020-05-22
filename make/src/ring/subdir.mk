
CPP_SRCS += \
./src/ring/Poly.cpp \
./src/ring/PolyRing.cpp 

OBJS += \
./make/src/ring/Poly.o \
./make/src/ring/PolyRing.o 

CPP_DEPS += \
./make/src/ring/Poly.d \
./make/src/ring/PolyRing.d 


make/src/ring/%.o: ./src/ring/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


