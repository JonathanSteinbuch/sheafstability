
CPP_SRCS += \
./src/jobs/JacobiTaylorJob.cpp \
./src/jobs/JobType.cpp \
./src/jobs/KernelJob.cpp \
./src/jobs/PowersJob.cpp \
./src/jobs/ReaderAndCaller.cpp \
./src/jobs/SemistabilityJob.cpp 

OBJS += \
./make/src/jobs/JacobiTaylorJob.o \
./make/src/jobs/JobType.o \
./make/src/jobs/KernelJob.o \
./make/src/jobs/PowersJob.o \
./make/src/jobs/ReaderAndCaller.o \
./make/src/jobs/SemistabilityJob.o 

CPP_DEPS += \
./make/src/jobs/JacobiTaylorJob.d \
./make/src/jobs/JobType.d \
./make/src/jobs/KernelJob.d \
./make/src/jobs/PowersJob.d \
./make/src/jobs/ReaderAndCaller.d \
./make/src/jobs/SemistabilityJob.d 


make/src/jobs/%.o: ./src/jobs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++17 -I/usr/local/include/ -O3 -fopenmp -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


