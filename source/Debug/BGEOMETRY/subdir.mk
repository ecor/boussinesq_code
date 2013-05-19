################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../BGEOMETRY/bigcells2.c \
../BGEOMETRY/g_raster2plvector.c \
../BGEOMETRY/geometry.c \
../BGEOMETRY/geometry2.c \
../BGEOMETRY/geometry_attribute.c \
../BGEOMETRY/geometry_freememory.c \
../BGEOMETRY/geometry_io.c \
../BGEOMETRY/geometry_utilities.c \
../BGEOMETRY/sorting.c 

OBJS += \
./BGEOMETRY/bigcells2.o \
./BGEOMETRY/g_raster2plvector.o \
./BGEOMETRY/geometry.o \
./BGEOMETRY/geometry2.o \
./BGEOMETRY/geometry_attribute.o \
./BGEOMETRY/geometry_freememory.o \
./BGEOMETRY/geometry_io.o \
./BGEOMETRY/geometry_utilities.o \
./BGEOMETRY/sorting.o 

C_DEPS += \
./BGEOMETRY/bigcells2.d \
./BGEOMETRY/g_raster2plvector.d \
./BGEOMETRY/geometry.d \
./BGEOMETRY/geometry2.d \
./BGEOMETRY/geometry_attribute.d \
./BGEOMETRY/geometry_freememory.d \
./BGEOMETRY/geometry_io.d \
./BGEOMETRY/geometry_utilities.d \
./BGEOMETRY/sorting.d 


# Each subdirectory must supply rules for building sources it contributes
BGEOMETRY/%.o: ../BGEOMETRY/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DWRITE_ITERATION_NUMBER -I../ASCII -I../BGEOMETRY -I../Boussinesq -I../FLUIDTURTLES -I../KeyPalette -I../MATH2 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


