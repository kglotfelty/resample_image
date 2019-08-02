tool = resample_image

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

resample_image_SOURCES = poly_clip.c resample_img.c t_resample_img.c
resample_image_CPPFLAGS = $(CIAO_CFLAGS)
resample_image_LDADD = $(CIAO_LIBS) -ldmimgio 
resample_image_LINK = $(CXX) -o $@ -Wl,-rpath,$(prefix)/lib -Wl,-rpath,$(prefix)/ots/lib 

if LINUX
resample_image_LDADD += -L$(prefix)/ots/lib -lstdc++
endif


dist_param_DATA = $(tool).par
#dist_ahelp_DATA = $(tool).xml

install-data-hook:
	chmod 444 $(paramdir)/$(dist_param_DATA)
#chmod 444 $(ahelpdir)/$(dist_ahelp_DATA)

