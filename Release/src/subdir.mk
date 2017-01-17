################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/RNAComplex.cpp \
../src/RNAGate_v18.cpp \
../src/RNAObject.cpp \
../src/Reaction.cpp \
../src/global_variables.cpp \
../src/nupack_conc_object.cpp \
../src/permutation_combination.cpp \
../src/random.cpp \
../src/read_write.cpp 

OBJS += \
./src/RNAComplex.o \
./src/RNAGate_v18.o \
./src/RNAObject.o \
./src/Reaction.o \
./src/global_variables.o \
./src/nupack_conc_object.o \
./src/permutation_combination.o \
./src/random.o \
./src/read_write.o 

CPP_DEPS += \
./src/RNAComplex.d \
./src/RNAGate_v18.d \
./src/RNAObject.d \
./src/Reaction.d \
./src/global_variables.d \
./src/nupack_conc_object.d \
./src/permutation_combination.d \
./src/random.d \
./src/read_write.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/usr/local/bin/g++ -fopenmp -m64 -I/usr/local/include -I/usr/local/include/c++/4.6.0 -I/usr/local/include/c++/4.6.0/backward -I/usr/local/lib/gcc/x86_64-apple-darwin10.6.0/4.6.0/include -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include" -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include/nupack" -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include/ViennaRNA/include/ViennaRNA" -O3 -ftree-vectorize -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


