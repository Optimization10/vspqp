#===============================================================================
# Demo/makefile
#===============================================================================


include $(VSPATH)/Config/internal.mk

.PHONY: example

example: $(VSPATH)/Demo/example


$(EXAMPLE): $(VSPATH)/Demo/example.c $(LIB_VSPQP)
	@echo ""
	@echo "Building executable $@..."
	$(CC) $(CFLAGS) $< $(LIB_VSPQP) $(LIBS) -o $@
	@echo "Success!"

$(LIB_VSPQP): FORCE
	$(MAKE) -C $(VSPATH)/Lib

FORCE:

-include $(VSPATH)/Demo/example.d