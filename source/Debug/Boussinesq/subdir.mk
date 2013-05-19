################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Boussinesq/b_readgrid.c \
../Boussinesq/b_solver.c \
../Boussinesq/b_sources.c \
../Boussinesq/b_utilities.c \
../Boussinesq/b_v_advection.c \
../Boussinesq/b_volumes.c \
../Boussinesq/main.c 

OBJS += \
./Boussinesq/b_readgrid.o \
./Boussinesq/b_solver.o \
./Boussinesq/b_sources.o \
./Boussinesq/b_utilities.o \
./Boussinesq/b_v_advection.o \
./Boussinesq/b_volumes.o \
./Boussinesq/main.o 

C_DEPS += \
./Boussinesq/b_readgrid.d \
./Boussinesq/b_solver.d \
./Boussinesq/b_sources.d \
./Boussinesq/b_utilities.d \
./Boussinesq/b_v_advection.d \
./Boussinesq/b_volumes.d \
./Boussinesq/main.d 


# Each subdirectory must supply rules for building sources it contributes
Boussinesq/%.o: ../Boussinesq/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


