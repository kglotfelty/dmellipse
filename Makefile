##############################################################################

###KJG MK_TOP = ../../../..
MK_TOP = /export/ciao_from_source/ciao-4.8/src
KJG = /export/ciao

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmellipse
LIB_FILES         =
PAR_FILES         = dmellipse.par
INC_FILES         =
XML_FILES         = dmellipse.xml

SRCS	= ellipse.c t_ellipse.c
OBJS	= $(SRCS:.c=.o)

LOCAL_LIBS = -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I$(MK_TOP)/da/analysis/dmtools/dmimgio/

MAKETEST_SCRIPT = dmellipse.t



include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LINK)
	@echo


kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |             Building dmellipse DM host tool           | "
	@echo "   \----------------------------------------------------------/ "
