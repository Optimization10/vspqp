#===============================================================================
# Demo/makefile
#===============================================================================


include $(VSPATH)/Config/internal.mk

.PHONY: example

example: $(VSPATH)/Demo/example


$(EXAMPLE): $(VSPATH)/Demo/example.c $(LIB_VSSTARTUP)
	@echo ""
	@echo "Building executable $@..."
	$(CC) $(CFLAGS) $< $(LIB_VSSTARTUP) $(LIBS) -o $@
	@echo "Success!"

$(LIB_VSSTARTUP): FORCE
	$(MAKE) -C $(VSPATH)/Lib

FORCE:

-include $(VSPATH)/Demo/example.d