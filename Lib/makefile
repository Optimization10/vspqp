#===============================================================================
# VS/Lib/Makefile: for compiling the static VSPQP library
#===============================================================================

include $(VSPATH)/Config/internal.mk

$(LIB_VSPQP): $(OBJS)
	@echo ""
	@echo "Building library $@..."
	ar rcs $@ $(OBJS)
	@echo "Success!"

-include $(OBJS:.o=.d)