################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ASCII/import_ascii.c \
../ASCII/rw_maps.c \
../ASCII/tabs.c \
../ASCII/write_ascii.c 

OBJS += \
./ASCII/import_ascii.o \
./ASCII/rw_maps.o \
./ASCII/tabs.o \
./ASCII/write_ascii.o 

C_DEPS += \
./ASCII/import_ascii.d \
./ASCII/rw_maps.d \
./ASCII/tabs.d \
./ASCII/write_ascii.d 


# Each subdirectory must supply rules for building sources it contributes
ASCII/%.o: ../ASCII/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


