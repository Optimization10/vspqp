# =============================================================================
# DO NOT EDIT - Internal build configuration
# Edit user_config.mk for user-configurable options
# =============================================================================

include $(VSPATH)/Config/user_config.mk


# =============================================================================
# Common Variables
# =============================================================================
SRC_DIR = $(VSPATH)/Source
BUILD_DIR = $(VSPATH)/Build

SRCS = $(shell find $(SRC_DIR) -name '*.c')
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

LIB_VSPQP = $(VSPATH)/Lib/libVSPQP.a
EXAMPLE = $(VSPATH)/Demo/example

# =============================================================================
# Basic Flags
# =============================================================================

ifeq ($(EXTRA_WARNINGS),enabled)
WARNINGS = -W -Wall -pedantic -Wmissing-prototypes \
	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
	-Wno-parentheses -Wshadow -Wcast-align \
	-Winline -Wstrict-prototypes -Wno-unknown-pragmas
endif

CFLAGS = -g3 $(WARNINGS) -I$(VSPATH)/Source/Include

# Used to handle dependency of object files on headers
CFLAGS += -MMD -MP

LIBS = -lm

# =============================================================================
# Debug or Release Version
# =============================================================================
ifeq ($(VERSION),debug)
CFLAGS += -O0 -DVS_DEBUG
else
CFLAGS += -O3 -flto -march=native
endif

# =============================================================================
# Types and MKL
# =============================================================================
ifeq ($(USE_MKL),true)
	CFLAGS += -DUSE_MKL -DVSFLOAT=double -I$(MKL_INCLUDE) -fopenmp
	LIBS += -lpthread -ldl

#   Use GNU threading layer for gcc, Intel threading layer for clang/icx
	ifneq (,$(findstring gcc,$(CC)))
		MKL_THREAD = -lmkl_gnu_thread
		MKL_OMP = -lgomp
	else
		MKL_THREAD = -lmkl_intel_thread
		MKL_OMP = -liomp5
	endif

	ifeq ($(LONG),true)
		CFLAGS += -DVSINT="long long int" -DVSLINT="long long int" -DMKL_ILP64
		MKL_LIBS = -Wl,--start-group -lmkl_intel_ilp64 $(MKL_THREAD) \
		-lmkl_core -Wl,--end-group $(MKL_OMP)
	else
		CFLAGS += -DVSINT=int -DVSLINT=int
		MKL_LIBS = -Wl,--start-group -lmkl_intel_lp64 $(MKL_THREAD) \
		-lmkl_core -Wl,--end-group $(MKL_OMP)
	endif

	LIBS := $(MKL_LIBS) $(LIBS)

else
#   Set the float type as defined by user, or to double if undefined
	ifndef VSFLOAT
		CFLAGS += -DVSFLOAT=double
	else
		CFLAGS += -DVSFLOAT=$(VSFLOAT)
	endif

#   Set the integer types
	ifeq ($(LONG),true)
		CFLAGS += -DVSINT="long long int" -DVSLINT="long long int"
	else
#       Define VSINT set by the user. If user did not set, default it to int
		ifndef VSINT
			CFLAGS += -DVSINT=int
		else
			CFLAGS += -DVSINT=$(VSINT)
		endif

#       Define VSLINT set by the user. If user did not set, default it to int
		ifndef VSLINT
			CFLAGS += -DVSLINT=int
		else
			CFLAGS += -DVSLINT=$(VSLINT)
		endif
	endif
endif

# =============================================================================
# Pattern Rules
# =============================================================================
# Build an object file from a source file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@echo ""
	@echo "Compiling $<..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Success!"