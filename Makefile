#Compilers
CXX			:=	g++
NVCC		:=	nvcc

# paths
INCLUDE_DIR		:= /global/scratch/users/seokjinmoon/organicSolvents/organicSolvents/include
SRC_DIR				:= src
LIB_DIR				:= /global/scratch/users/seokjinmoon/organicSolvents/organicSolvents/src
BUILD_DIR   	:= /global/scratch/users/seokjinmoon/organicSolvents/organicSolvents/build
BIN_DIR 			:= bin

# Flags
CXXFLAGS  		:= -std=c++17 -O2 -I$(INCLUDE_DIR)
NVCCFLAGS 		:= -std=c++17 -O2 -I$(INCLUDE_DIR)

# ------ Library sources ------
LIB_SRCS			:= $(wildcard $(LIB_DIR)/*.cpp)
LIB_OBJS			:= $(patsubst $(LIB_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(LIB_SRCS))

#programs to build
PROGRAMS := $(notdir $(wildcard $(SRC_DIR)/*))

# Helper
define SRC_FILES
$(wildcard $(SRC_DIR)/$1/*.cu)
endef

# Target executable
TARGETS = $(patsubst %, $(BIN_DIR)/%, $(PROGRAMS))

.PHONY: all clean
all: $(TARGETS)

# Ensure build/output directories exists
$(BUILD_DIR):
	@mkdir -p $@
$(BIN_DIR):
	@mkdir -p $@

# Build library
$(BUILD_DIR)/%.o: $(LIB_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

## create map: executable depends on its source files
define PROGRAM_RULE
$(BIN_DIR)/$1: $$(call SRC_FILES,$1) $(LIB_OBJS) | $(BIN_DIR)
	@echo "Building $1 ..."
	$$(NVCC) $$(NVCCFLAGS) $$(call SRC_FILES,$1) $$(LIB_OBJS) -o $$@
endef


# one build rule per program
$(foreach prog,$(PROGRAMS), $(eval $(call PROGRAM_RULE,$(prog))))

# allow "make program1" instead of "make bin/program1"
$(PROGRAMS): %: $(BIN_DIR)/%
	@echo "Done building $@"

clean:
	rm -rf $(BIN_DIR)


