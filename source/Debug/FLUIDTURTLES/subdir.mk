################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../FLUIDTURTLES/alloc.c \
../FLUIDTURTLES/datamanipulation.c \
../FLUIDTURTLES/error.c \
../FLUIDTURTLES/linearalgebra.c \
../FLUIDTURTLES/list.c \
../FLUIDTURTLES/probability.c \
../FLUIDTURTLES/random.c \
../FLUIDTURTLES/statistics.c \
../FLUIDTURTLES/string.c \
../FLUIDTURTLES/t_io.c \
../FLUIDTURTLES/tensors3D.c \
../FLUIDTURTLES/utilities.c \
../FLUIDTURTLES/write_dem.c 

OBJS += \
./FLUIDTURTLES/alloc.o \
./FLUIDTURTLES/datamanipulation.o \
./FLUIDTURTLES/error.o \
./FLUIDTURTLES/linearalgebra.o \
./FLUIDTURTLES/list.o \
./FLUIDTURTLES/probability.o \
./FLUIDTURTLES/random.o \
./FLUIDTURTLES/statistics.o \
./FLUIDTURTLES/string.o \
./FLUIDTURTLES/t_io.o \
./FLUIDTURTLES/tensors3D.o \
./FLUIDTURTLES/utilities.o \
./FLUIDTURTLES/write_dem.o 

C_DEPS += \
./FLUIDTURTLES/alloc.d \
./FLUIDTURTLES/datamanipulation.d \
./FLUIDTURTLES/error.d \
./FLUIDTURTLES/linearalgebra.d \
./FLUIDTURTLES/list.d \
./FLUIDTURTLES/probability.d \
./FLUIDTURTLES/random.d \
./FLUIDTURTLES/statistics.d \
./FLUIDTURTLES/string.d \
./FLUIDTURTLES/t_io.d \
./FLUIDTURTLES/tensors3D.d \
./FLUIDTURTLES/utilities.d \
./FLUIDTURTLES/write_dem.d 


# Each subdirectory must supply rules for building sources it contributes
FLUIDTURTLES/%.o: ../FLUIDTURTLES/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


