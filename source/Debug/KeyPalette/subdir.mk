################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../KeyPalette/additional_read_functions.c \
../KeyPalette/get_filenames.c \
../KeyPalette/key.palette.c \
../KeyPalette/read_command_line.c 

OBJS += \
./KeyPalette/additional_read_functions.o \
./KeyPalette/get_filenames.o \
./KeyPalette/key.palette.o \
./KeyPalette/read_command_line.o 

C_DEPS += \
./KeyPalette/additional_read_functions.d \
./KeyPalette/get_filenames.d \
./KeyPalette/key.palette.d \
./KeyPalette/read_command_line.d 


# Each subdirectory must supply rules for building sources it contributes
KeyPalette/%.o: ../KeyPalette/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


