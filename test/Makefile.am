TESTS = resample_image.t

TESTS_ENVIRONMENT =\
TESTIN=indata; \
TESTSAV=save_data; \
TESTOUT=delme; \
TESTLOG=delme; \
TESTRUN=delme; \
PFILES=$(TESTOUT)\;$(top_srcdir)/src:$(prefix)/param:$(prefix)/contrib/param; \
PATH=$(top_srcdir)/src:$(prefix)/bin:$(prefix)/ots/bin:$(PATH); \
ASCDS_INSTALL=$(prefix); \
export TESTIN TESTSAV TESTOUT TESTLOG TESTRUN PFILES LD_LIBRARY_PATH ASCDS_INSTALL; 

clean-local:
	-rm -rf delme
