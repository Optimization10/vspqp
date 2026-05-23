#===============================================================================
# Config/user_config.mk
#===============================================================================

# NOTE: Run "make clean" after changing any configurations in this file

#------------------------------------------------------------------------------
# Choose debug or release version
#------------------------------------------------------------------------------
VERSION = release
# VERSION = dev
# VERSION = debug


#------------------------------------------------------------------------------
# Configure Intel MKL
#------------------------------------------------------------------------------

# Setting USE_MKL = true enables Intel MKL
# Comment out the line below to disable Intel MKL
# USE_MKL = true

# Path to Intel MKL include files. Ignored if MKL not used
MKL_INCLUDE = /usr/include/mkl
 
#------------------------------------------------------------------------------
# Configure Types
#------------------------------------------------------------------------------

# If MKL is enabled, this choice is ignored and VSPFLOAT is set to double
VSPFLOAT = double
# VSPFLOAT = long double


# If LONG = true, both VSPINT and VSPLONG are set to 64-bit integers.
# If LONG is not set:
#   - With MKL enabled, VSPINT and VSPLONG default to 32-bit integers
#   - With MKL disabled, VSPINT and VSPLONG can be set independently below


# ---- Used only if MKL is disabled and LONG is not true ----------------------
VSPINT = int
# VSPINT = long long int

VSPLONG = int
# VSPLONG = long long int

#------------------------------------------------------------------------------
# Enable or disable extra compiler warnings
#------------------------------------------------------------------------------
# EXTRA_WARNINGS = enabled
EXTRA_WARNINGS = disabled

#------------------------------------------------------------------------------
# Specify C Compiler
#------------------------------------------------------------------------------
CC = gcc
# CC = clang


#------------------------------------------------------------------------------
# Additional Compiler Flags
#------------------------------------------------------------------------------
# Add any other compiler flags, such as as -DLOG
# USER_FLAGS = -fsanitize=address