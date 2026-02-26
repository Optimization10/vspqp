# --- Complete Makefile for VSPQP Project ---

# 1. Configuration
# ------------------------------------
TARGET = build/vspqp
CC = gcc
LDFLAGS = -lm -lrt -L$(HOME)/.local/lib -Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lpthread -ldl -lm
BUILD_DIR = build

# --- INTEGER MODE CONFIGURATION ---
# Default mode for VSP_INTEGERS_MODE is 0 if not explicitly set
# 0 --> VSPINT = int        and  VSPLONG = int
# 1 --> VSPINT = int        and  VSPLONG = long long
# 2 --> VSPINT = long long  and  VSPLONG = long long
VSP_INTEGERS_MODE = 0

# Source and Include Directories (ensure these paths are correct)
SRCDIRS   = source
INCDIRS   = include source $(HOME)/.local/include /usr/include/mkl
OUTPUTDIR = output

# Define common flags and optimization/debug specific flags
COMMON_CFLAGS = -g -DVSP_INTEGERS_MODE=$(VSP_INTEGERS_MODE)
DEFINES = -DMEASURE_PERFORMANCE

# Performance flags: Optimization, native architecture, Link-Time Optimization (LTO)
CFLAGS_PERFORMANCE = $(COMMON_CFLAGS) -O3 -march=native -flto $(DEFINES)

# Debug flags: No optimization, full debugging symbols
CFLAGS_DEBUG = $(COMMON_CFLAGS) -O0 -g3 -DVSP_DEBUG

# VS code flags: no optimizations, but no -DVSP_DEBUG either
CFLAGS_VSCODE = $(COMMON_CFLAGS) -O0 -g3 

# Valgrind flags
VALGRIND_FLAGS = --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose

# Address Sanitizer flags: No optimization (-O0), full symbols (-g3), plus sanitizer flag
CFLAGS_ASAN = $(COMMON_CFLAGS) -O0 -g3 -fsanitize=address,undefined,bounds


# 2. Recursive File and Directory Discovery
# ------------------------------------
# Find all .c files in SRCDIRS and their subdirectories
C_SOURCES := $(shell find $(SRCDIRS) -name "*.c")

# Find ALL unique directories that contain *.h files, recursively.
# This solves the non-path-qualified includes issue by passing every header directory
# to the compiler with -I.
HEADER_DIRS := $(shell find $(SRCDIRS) $(INCDIRS) -type f -name "*.h" | xargs -n1 dirname | sort -u)

# Create the -I flags from the discovered directories
INCLUDES = $(addprefix -I,$(HEADER_DIRS))

# Convert source file list (.c) to object file list (.o) in the build directory
OBJECTS := $(patsubst %.c,$(BUILD_DIR)/%.o,$(C_SOURCES))
DEPFILES := $(OBJECTS:.o=.d)


# 3. Compiler Setup
# ------------------------------------
# Define the DEFAULT CFLAGS to be used if no specific target is called (e.g., 'make' alone)
# The '?=' operator assigns only if CFLAGS is not already set.
CFLAGS ?= $(CFLAGS_DEBUG)

# Append the -I flags to the chosen CFLAGS for the current build mode
CFLAGS += $(INCLUDES)


# 4. Targets
# ------------------------------------
.PHONY: all clean performance debug valgrind asan

# Default target: build the debug version
all: debug

# Performance Target
# This target temporarily overrides CFLAGS to use the performance flags
performance: CFLAGS = $(CFLAGS_PERFORMANCE) $(INCLUDES)
performance: $(TARGET)
	@echo "Performance Build Complete."
	@./$(TARGET)

# Debug Target
# This target temporarily overrides CFLAGS to use the debug flags
debug: CFLAGS = $(CFLAGS_DEBUG) $(INCLUDES)
debug: $(TARGET)
	@echo "Debug Build Complete."
	./$(TARGET)

# VS Code Target
# This target temporarily overrides CFLAGS to use the debug flags
vscode_dbg: CFLAGS = $(CFLAGS_VSCODE) $(INCLUDES) -DVSP_DEBUG
vscode_dbg: $(TARGET)
	@echo "Debug Build Complete."

# Valgrind Target: Builds with NO VSP_DEBUG, then runs through Valgrind
valgrind: CFLAGS = $(CFLAGS_VALGRIND) $(INCLUDES)
valgrind: $(TARGET)
	@mkdir -p output # Ensure the directory exists
	@echo "Running Valgrind Memcheck on $(TARGET). Output redirected to output/valgrind.txt"
	@valgrind $(VALGRIND_FLAGS) ./$(TARGET) 2> output/valgrind.txt

# Address Sanitizer Target
# Builds with multiple sanitizers and redirects output to output/asan.txt
asan: CFLAGS = $(CFLAGS_ASAN) $(INCLUDES)
asan: $(TARGET)
	@mkdir -p output # Ensure the directory exists
	@echo "Running ASan checks on $(TARGET). Output redirected to output/asan.txt"
	@ASAN_OPTIONS=abort_on_error=1 ./$(TARGET) 2> output/asan.txt

# Link the executable
# Note: This is a secondary prerequisite for the performance/debug targets.
$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(CC) $(OBJECTS) $(CFLAGS) -o $@ $(LDFLAGS)


# Include generated dependency files (critical for tracking header changes)
# The '-' prevents make from failing if the files don't exist yet (first run)
-include $(DEPFILES)

# Compile Sources (Pattern Rule)
# This rule is used for every object file and uses the currently active CFLAGS
$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(BUILD_DIR)/$(dir $<)
	@$(CC) $(CFLAGS) -c $< -o $@ -MMD -MP

# Clean up the build files
clean:
	@rm -rf $(BUILD_DIR) $(TARGET)
	@rm -f $(shell find . -name '*.d')
	@rm ./output/*
	@echo "Cleaned build artifacts."