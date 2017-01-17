################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/RNAComplex.cpp \
../src/RNAObject.cpp \
../src/Reaction.cpp \
../src/global_variables.cpp \
../src/main.cpp \
../src/nupack_conc_object.cpp \
../src/permutation_combination.cpp \
../src/random.cpp \
../src/read_write.cpp 

OBJS += \
./src/RNAComplex.o \
./src/RNAObject.o \
./src/Reaction.o \
./src/global_variables.o \
./src/main.o \
./src/nupack_conc_object.o \
./src/permutation_combination.o \
./src/random.o \
./src/read_write.o 

CPP_DEPS += \
./src/RNAComplex.d \
./src/RNAObject.d \
./src/Reaction.d \
./src/global_variables.d \
./src/main.d \
./src/nupack_conc_object.d \
./src/permutation_combination.d \
./src/random.d \
./src/read_write.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	/usr/local/bin/g++ -I/usr/local/include -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include" -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include/nupack" -I"/Users/kamal/Documents/my_programs/c_code/RNAGate_v18/include/ViennaRNA/include/ViennaRNA" -O2 -g -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


