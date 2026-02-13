# Compiler and flags
ifeq ($(CPU), 1)
    CXX = g++
    CXXFLAGS += -std=c++17 -O3
    ifeq ($(DEBUG), 1)
        CXXFLAGS += -g
    endif
else
    CXX = nvcc
    CXXFLAGS += -O3 -std=c++17 -DUSE_CUDA
    ifeq ($(DEBUG), 1)
        CXXFLAGS += -G -g
    endif
endif

CXXFLAGS += -Iinclude -Isrc/precip -Isrc/h2spec -Isrc/world1d -Isrc/utils -Isrc/iondens -Isrc/h3spec

# Output dirs
OBJDIR := obj
BINDIR := bin

# Source files
CPP_SRCS := $(shell find src -name "*.cpp")
CU_SRCS  := $(shell find src -name "*.cu")

# Object files (mirror src/ path under obj/)
CPP_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPP_SRCS))
CU_OBJS  := $(patsubst %.cu,$(OBJDIR)/%.o,$(CU_SRCS))

ifeq ($(CPU), 1)
    OBJS := $(CPP_OBJS)
else
    OBJS := $(CPP_OBJS) $(CU_OBJS)
endif

# Dependency files (for .cpp only; nvcc dep-gen varies by setup)
DEPS := $(CPP_OBJS:.o=.d)

# Executable name/location
EXEC := $(BINDIR)/jaic

# Default target
all: $(EXEC)

# Ensure output directories exist
$(BINDIR) $(OBJDIR):
	mkdir -p $@

# Compile .cpp -> obj/...
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Compile .cu -> obj/...
$(OBJDIR)/%.o: %.cu | $(OBJDIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link
$(EXEC): $(OBJS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

# Clean up
clean:
	rm -rf $(OBJDIR) $(BINDIR)

# Include dependency files (if they exist)
-include $(DEPS)

.PHONY: all clean
