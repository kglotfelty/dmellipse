##############################################################################

###KJG MK_TOP = ../../../..
MK_TOP = /export/ciao_from_source/ciao-4.7/src

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmellipse
LIB_FILES         =
PAR_FILES         = dmellipse.par
INC_FILES         =
XML_FILES         = dmellipse.xml

SRCS	= ellipse.c t_ellipse.c
OBJS	= $(SRCS:.c=.o)

###KJG LOCAL_LIBS = -L../dmimgio/ -ldmimgio
###KJG LOCAL_INC  = -I../dmimgio/

LOCAL_LIBS = -L/export/ciao_from_source/ciao-4.7/src/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I/export/ciao_from_source/ciao-4.7/src/da/analysis/dmtools/dmimgio/

MAKETEST_SCRIPT = dmellipse.t



include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LINK)
	@echo

announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |             Building dmellipse DM host tool           | "
	@echo "   \----------------------------------------------------------/ "
