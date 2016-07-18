
#MK_TOP = ../../../..
MK_TOP = /export/ciao_from_source/ciao-4.8/src
KJG = /export/ciao

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = resample_image
PURE_EXEC	  = 
LIB_FILES         =
PAR_FILES         = resample_image.par
INC_FILES         =
XML_FILES         = 

SRCS	= resample_img.c t_resample_img.c poly_clip.c

OBJS	= $(SRCS:.c=.o)


MAKETEST_SCRIPT = r

LOCAL_INC 	= -I$(MK_TOP)/da/analysis/dmtools/dmimgio/
LOCAL_LIBS 	= -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio



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
	@echo "   |             Building resample_image DM host tool         | "
	@echo "   \----------------------------------------------------------/ "
