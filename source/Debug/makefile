################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include MATH2/subdir.mk
-include KeyPalette/subdir.mk
-include FLUIDTURTLES/subdir.mk
-include Boussinesq/subdir.mk
-include BGEOMETRY/subdir.mk
-include ASCII/subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: 0.9375KMackenzie_Boussinesq_svn

# Tool invocations
0.9375KMackenzie_Boussinesq_svn: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C Linker'
	gcc  -o "0.9375KMackenzie_Boussinesq_svn" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) 0.9375KMackenzie_Boussinesq_svn
	-@echo ' '

post-build:
	-cp 0.9375KMackenzie_Boussinesq_svn /Users/ecor/local/bin/boussinesq
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: post-build

-include ../makefile.targets
