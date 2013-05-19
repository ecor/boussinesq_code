################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../MATH2/boussinesq_matrix.c \
../MATH2/doublevector_utilities.c \
../MATH2/geo_statistic.09375.c \
../MATH2/linear_span.c \
../MATH2/preconditioned_conjugate_gradient.c \
../MATH2/util_math.c 

OBJS += \
./MATH2/boussinesq_matrix.o \
./MATH2/doublevector_utilities.o \
./MATH2/geo_statistic.09375.o \
./MATH2/linear_span.o \
./MATH2/preconditioned_conjugate_gradient.o \
./MATH2/util_math.o 

C_DEPS += \
./MATH2/boussinesq_matrix.d \
./MATH2/doublevector_utilities.d \
./MATH2/geo_statistic.09375.d \
./MATH2/linear_span.d \
./MATH2/preconditioned_conjugate_gradient.d \
./MATH2/util_math.d 


# Each subdirectory must supply rules for building sources it contributes
MATH2/%.o: ../MATH2/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


