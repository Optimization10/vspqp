# NOTE: Run "make clean" after changing any configurations in this file

#------------------------------------------------------------------------------
# Choose debug or release version
#------------------------------------------------------------------------------
VERSION = debug
# VERSION = release

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

# If MKL is enabled, this choice is ignored and VSFLOAT is set to double
VSFLOAT = double
# VSFLOAT = long double


# If LONG = true, both VSINT and VSLINT are set to 64-bit integers.
# If LONG is not set:
#   - With MKL enabled, VSINT and VSLINT default to 32-bit integers
#   - With MKL disabled, VSINT and VSLINT can be set independently below


# ---- Used only if MKL is disabled and LONG is not true ----------------------
VSINT = int
# VSINT = long long int

VSLINT = int
# VSLINT = long long int

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
