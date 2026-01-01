
CXX = clang++
LOCAL_INCLUDE = /Users/xtao/local/include
HDF5_INCLUDE = /opt/local/include
HDF5_LIB = /opt/local/lib

OPT = -O2

SRC_DIR := source
BUILD_DIR = build

DIRS := $(shell find $(SRC_DIR) -type d)
SRCS := $(shell find $(SRC_DIR)/* -name \*.cc)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cc=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cc=.d))

# OpenMP (MacPorts example; adjust if needed)
OPENMP_CFLAGS = -Xpreprocessor -fopenmp -I/opt/local/include/libomp
OPENMP_LDFLAGS = -L/opt/local/lib/libomp
OPENMP_LIBS = -lomp

CXXFLAGS = $(OPT) -I$(HDF5_INCLUDE) -I$(LOCAL_INCLUDE) $(DIRS:%=-I%) $(OPENMP_CFLAGS)
LDFLAGS  = -L$(HDF5_LIB) $(OPENMP_LDFLAGS)
LDLIBS   = -lhdf5 $(OPENMP_LIBS)

executable = sayram-2d

.PHONY: all clean
all: $(executable)

$(executable): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

$(BUILD_DIR)/%.d: %.cc
	@echo "Checking dependencies for $<"
	@if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(CXX) $(CXXFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

$(BUILD_DIR)/%.o: %.cc
	@echo "Compiling $<"
	@if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "Cleaning $(BUILD_DIR)"
	rm -r $(BUILD_DIR)
