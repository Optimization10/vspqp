# =============================================================================
# Config/internal.mk
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
CFLAGS = $(USER_FLAGS)

ifeq ($(EXTRA_WARNINGS),enabled)
WARNINGS = -W -Wall -pedantic -Wmissing-prototypes \
	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
	-Wno-parentheses -Wshadow -Wcast-align \
	-Winline -Wstrict-prototypes -Wno-unknown-pragmas
endif

CFLAGS += -g3 $(WARNINGS) -I$(VSPATH)/Include -I$(VSPATH)/Source/Include

# Used to handle dependency of object files on headers
CFLAGS += -MMD -MP

LIBS = -lm

# =============================================================================
# Debug or Release Version
# =============================================================================
ifeq ($(VERSION),debug)
CFLAGS += -O0 -DDEBUG
else ifeq ($(VERSION),dev)
CFLAGS += -O0
else
CFLAGS += -O3 -flto -march=native
endif

# =============================================================================
# Types and MKL
# =============================================================================
ifeq ($(USE_MKL),true)
	CFLAGS += -DUSE_MKL -DVSPFLOAT=double -I$(MKL_INCLUDE) -fopenmp
	RESOLVED_VSPFLOAT=double
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
		CFLAGS += -DVSPINT="long long int" -DVSPLONG="long long int" -DMKL_ILP64
		RESOLVED_VSPINT=long long int
		RESOLVED_VSPLONG=long long int
		MKL_LIBS = -Wl,--start-group -lmkl_intel_ilp64 $(MKL_THREAD) \
		-lmkl_core -Wl,--end-group $(MKL_OMP)
	else
		CFLAGS += -DVSPINT=int -DVSPLONG=int
		RESOLVED_VSPINT=int
		RESOLVED_VSPLONG=int
		MKL_LIBS = -Wl,--start-group -lmkl_intel_lp64 $(MKL_THREAD) \
		-lmkl_core -Wl,--end-group $(MKL_OMP)
	endif

	LIBS := $(MKL_LIBS) $(LIBS)

else
#   Set the float type as defined by user, or to double if undefined
	ifndef VSPFLOAT
		CFLAGS += -DVSPFLOAT=double
		RESOLVED_VSPFLOAT=double
	else
		CFLAGS += -DVSPFLOAT=$(VSPFLOAT)
		RESOLVED_VSPFLOAT=$(VSPFLOAT)
	endif

#   Set the integer types
	ifeq ($(LONG),true)
		CFLAGS += -DVSPINT="long long int" -DVSPLONG="long long int"
		RESOLVED_VSPINT=long long int
		RESOLVED_VSPLONG=long long int
	else
#       Define VSPINT set by the user. If user did not set, default it to int
		ifndef VSPINT
			CFLAGS += -DVSPINT=int
			RESOLVED_VSPINT=int
		else
			CFLAGS += -DVSPINT=$(VSPINT)
			RESOLVED_VSPINT=$(VSPINT)
		endif

#       Define VSPLONG set by the user. If user did not set, default it to int
		ifndef VSPLONG
			CFLAGS += -DVSPLONG=int
			RESOLVED_VSPLONG=int
		else
			CFLAGS += -DVSPLONG=$(VSPLONG)
			RESOLVED_VSPLONG=$(VSPLONG)
		endif
	endif
endif

# =============================================================================
# Confid Header Path
# =============================================================================
VSPQP_CONFIG_H = $(VSPATH)/Include/VSPQP_config.h


# =============================================================================
# Pattern Rules
# =============================================================================
# Build an object file from a source file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(VSPQP_CONFIG_H)
	@echo ""
	@echo "Compiling $<..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Success!"
