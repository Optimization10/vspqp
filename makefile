#===============================================================================
# Main Makefile
#===============================================================================
VSPATH = $(realpath $(CURDIR))
export VSPATH

include $(VSPATH)/Config/internal.mk

.PHONY: lib example run_example valgrind clean

default: lib

lib: 
	$(MAKE) -C $(VSPATH)/Lib


# To make a new example (such as Demo/X.c)
#    In the current makefile:
#    	1) Add X target to this file similar to the example target
#       2) and add X to .PHONY above
#    	3) Add rm -f lines for X and X.d in clean target
#       4) Add run_X target similar to run_example
#    
#	 In the demo makefile:
#	    1) Add X to .PHONY
#       2) Add X: $(VSPATH)/Demo/X
#       3) Add -include $(VSPATH)/Demo/X.d
example:
	$(MAKE) -C $(VSPATH)/Demo example





run_example: example
	@echo ""
	@echo "Running example...\n"
	@./Demo/example

valgrind: example
	@valgrind \
  --tool=memcheck \
  --leak-check=full \
  --show-leak-kinds=all \
  --track-origins=yes \
  --suppressions=$(VSPATH)/Config/mkl.supp \
  ./Demo/example

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(LIB_VSPQP)
	rm -f $(VSPATH)/Demo/example
	rm -f $(VSPATH)/Demo/example.d
