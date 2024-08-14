# Compiler information; this makefile is based on the one for Smilie project
# This Makefile used be used in the parent folder of source ("../")

CC = g++
LOCAL_INCLUDE = /Users/xtao/local/include # change to your own include path
# HDF5_LIB = /Users/username/local/lib
SRC_DIR := source
BUILD_DIR = build

DIRS := $(shell find $(SRC_DIR) -type d)
SRCS := $(shell find $(SRC_DIR)/* -name \*.cc)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cc=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cc=.d))

CCFLAGS = -Wall -Wno-class-memaccess -O2 -fopenmp -I$(LOCAL_INCLUDE) 
CCFLAGS += $(DIRS:%=-I%)
CCFLAGS += 

# LDFLAGS = -L$(HDF5_LIB) -lhdf5
LDFLAGS = 

executable= fvm2d

.PHONY: all clean

#-----------------------------------------------------
# Set the verbosity prefix
ifeq (,$(findstring verbose,$(config)))
    Q := @
else
    Q :=
endif

all:  $(executable)

# link objs
$(executable):$(OBJS) 
	$(CC) $(LDFLAGS) $(OBJS) -o $@

# dependences 
$(BUILD_DIR)/%.d: %.cc
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(CC) $(CCFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

# objects 
$(BUILD_DIR)/%.o: %.cc
	@echo "Compiling $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(CC) $(CCFLAGS) -c $< -o $@

clean:
	@echo "Cleaning $(BUILD_DIR)"
	$(Q) rm -r $(BUILD_DIR)

