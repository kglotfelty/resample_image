TESTS = resample_image.t

TESTS_ENVIRONMENT =\
TESTIN=indata; \
TESTSAV=save_data; \
TESTOUT=delme; \
TESTLOG=delme; \
TESTRUN=delme; \
export TESTIN TESTSAV TESTOUT TESTLOG TESTRUN; \
PATH=$(PATH):$(top_srcdir)/src;\
PFILES=$(TESTOUT)\;$(top_srcdir)/src:$(prefix)/param:$(prefix)/contrib/param; \
LD_LIBRARY_PATH=$(prefix)/lib:$(prefix)/ots/lib \
source $(prefix)/bin/ciao.sh;

clean-local:
	-rm -rf delme
