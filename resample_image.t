#! /bin/sh

# 30 January 2002

# This is the official template for pipetool regression test scripts.
# In addition to supporting the "SHORTTEST" option, this script also
# allows the user to run individual subtests from the command line.
# The script will accept a series of test identifiers, generally of
# the form "test1" "test2" ... which are to be run.

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!3, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!4 below.


# !!1
# resample_image.t
# test script for resample_image


# !!2
# syntax:
# resample_image.t [<testid> ... ]
 



######################################################################
# subroutine
# error_exit <message>
# Fatal error exit

error_exit()
{
  echo "$1" | tee -a $LOGFILE
  echo "${toolname} : FAIL" | tee -a $LOGFILE
  exit 1
}

######################################################################
# subroutine
# keyfilter infile outfile
# filters out CHECKSUM, Dataset, CREATOR, HISTORY, DATASUM, 
#             ASCDSVER, HISTNUM, and DATE
# To filter additional keywords, add s/KEYWORD/Dataset/g; for each.

keyfilter()
{
  cat $1 | sed -e 's/CHECKSUM/Dataset/g;s/COMMENT/Dataset/g;
  s/DATE/Dataset/g;s/CREATOR/Dataset/g;s/HISTORY/Dataset/g;
  s/DATASUM/Dataset/g;s/ASCDSVER/Dataset/g;s/HISTNUM/Dataset/g' | \
  grep -v Dataset > $2
  zerotest $2
}

######################################################################
# subroutine
# find_tool <toolname>
# checks that tool exists and is runnable

find_tool()
{
  s1=`type $1`
  s2=`echo $s1 | awk -F" " '{ print $3}'`
  if test -x $s2 ; then
    :
  else
    error_exit "tool $1 not found"
  fi
}

######################################################################
# subroutine
# zerotest <file> 
# Makes sure that file is not 0 length.
# Use this to protect yourself against empty files  (which will 
# 'diff' without error).  This can happen when the input file to
# cat $infile | do_something >> $outfile
# is missing.  This is used by keyfilter(), above.

zerotest()
{
 if test -s $1 ;
 then
   :
 else
   echo "ERROR: file $1 is of zero length" >> $LOGFILE
   #  Indicate failure, but do not exit.
   mismatch=0
 fi
}


######################################################################
# Initialization

# !!3
toolname="resample_image"

# set up list of tests
# !!4
alltests="test_full test_quick test_physical test_logical test_average full_calcgti"
# test_hst test_hst2 test_2mass"

# "short" test to run
# !!5
shortlist="$alltests"


# compute date string for log file
DT=`date +'%d%b%Y_%T'`


# convenience definitions
OUTDIR=$TESTOUT/$toolname
SAVDIR=$TESTSAV/$toolname
INDIR=$TESTIN/$toolname
LOGDIR=$TESTLOG/$toolname

# set up log file name
LOGFILE=$LOGDIR/${toolname}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${toolname}_log.*



# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi


# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR 
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined" 
fi


# check for tools
# if a utility is used in the form "utility <args> > outfile", and 'utility'
# cannot be run, 'outfile' will still be created.  If utility is used on 
# both the output and reference files of a tool the resultant utility output 
# files will both exist and be empty, and will pass a diff.

find_tool dmdiff
find_tool dmimgcalc



# announce ourselves
echo ""
echo "${toolname} regression" | tee $LOGFILE
echo ""

# All parameters except verbose should be set anyway, but clear them
# to be safe.
bose=`pget $toolname verbose`
punlearn $toolname
pset $toolname verbose=$bose

lkTab1=${ASCDS_CALIB}/dmmerge_header_lookup.txt
lkTab2=${INDIR}/calcgti_dmmerge_header_lookup.txt

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do
    
  # delete old outputs
  rm -f $OUTDIR/${testid}*

  # Set up file names
  outfile=$OUTDIR/${testid}.fits
  savfile=$SAVDIR/${testid}.fits

  echo "running $testid" >> $LOGFILE

  ####################################################################
  # run the tool
  case ${testid} in
    # !!6
    test_full ) test1_string="resample_image infile=$INDIR/acisf00934N002_full_img2.fits matchfile=$INDIR/acisf04731N001_full_img2.fits outfile=$outfile res=1 clob+ coord=world lookupTab=${lkTab1}"
            ;;

    # same as test_full except using calcgti/lkTab2
    full_calcgti ) test1_string="resample_image infile=$INDIR/acisf00934N002_full_img2.fits matchfile=$INDIR/acisf04731N001_full_img2.fits outfile=$outfile res=1 clob+ coord=world lookupTab=${lkTab2}"
            ;;

    test_quick ) test1_string="resample_image infile=$INDIR/acisf00934N002_full_img2.fits matchfile=$INDIR/acisf04731N001_full_img2.fits outfile=$outfile res=0 clob+ coord=world lookupTab=${lkTab1}"
            ;;

    test_physical ) test1_string="resample_image infile=$INDIR/acisf04731N001_full_img2.fits matchfile=$INDIR/bin16.fits outfile=$outfile res=1 clob+ coord=physical lookupTab=${lkTab1}"
            ;;

    test_logical ) test1_string="resample_image infile=$INDIR/acisf00934N002_full_img2.fits matchfile=$INDIR/acisf04731N001_full_img2.fits outfile=$outfile res=1 clob+ coord=logical lookupTab=${lkTab1}"
            ;;

    test_average ) test1_string="resample_image infile=$INDIR/acisf00934N002_full_img2.fits matchfile=$INDIR/acisf04731N001_full_img2.fits outfile=$outfile res=0 clob+ coord=world meth=average lookupTab=${lkTab1}"
            ;;

    test_hst ) test1_string="resample_image infile=$INDIR/orion_hst.fits matchfile=$INDIR/orion_chandra.fits outfile=$outfile res=1 clob+ coord=world meth=sum lookupTab=${lkTab1}"
            ;;

#
# same as above but use hst as match file, chandra as ref.
#

    test_hst2 ) test1_string="resample_image matchfile=$INDIR/orion_hst.fits infile=$INDIR/orion_chandra.fits outfile=$outfile res=1 clob+ coord=world meth=sum lookupTab=${lkTab1}"
            ;;

    test_2mass ) test1_string="resample_image infile=$INDIR/orion_2mass.fits matchfile=$INDIR/orion_chandra.fits outfile=$outfile res=1 clob+ coord=world meth=sum lookupTab=${lkTab1}"
            ;;

 
  esac

  echo "ciaorun ./"$test1_string | tee -a  $LOGFILE 
  eval "ciaorun ./"$test1_string  | tee -a  $LOGFILE  2>&1
 

  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1

  # if different tests need different kinds of comparisons, use a 
  #  case ${testid} in...  here

  ####################################################################
  # FITS table    (duplicate for as many tables per test as needed)

  # new output
  # !!10
#dmlist $outfile header,data,array > $OUTDIR/${testid}.dmp1  2>>$LOGFILE
#keyfilter $OUTDIR/${testid}.dmp1 $OUTDIR/${testid}.dmp2  2>>$LOGFILE

  # reference output
  # !!11
#dmlist $savfile header,data,array > $OUTDIR/${testid}.dmp1_std  2>>$LOGFILE
#keyfilter $OUTDIR/${testid}.dmp1_std $OUTDIR/${testid}.dmp2_std \
#          2>>$LOGFILE

  # compare
  # !!12
#diff $OUTDIR/${testid}.dmp2 $OUTDIR/${testid}.dmp2_std > \
#     /dev/null 2>>$LOGFILE


dmdiff $outfile $savfile tol=$SAVDIR/tolerance > /dev/null 2>>$LOGFILE
if  test $? -ne 0 ; then
  echo "ERROR: MISMATCH in $outfile" >> $LOGFILE
  mismatch=0
fi


  ####################################################################
  # FITS image  (duplicate for as many images per test as needed)

  # check image
  # !!13
#   dmimgcalc "$outfile[1]" "$savfile[1]" none tst verbose=0   2>>$LOGFILE
#   if test $? -ne 0; then
#     echo "ERROR: DATA MISMATCH in $outfile" >> $LOGFILE
#     mismatch=0
#   fi

  #  Check the header of the image

  # !!14
  # dmlist $outfile header > $OUTDIR/${testid}.dmp1  2>>$LOGFILE
  # keyfilter $OUTDIR/${testid}.dmp1 $OUTDIR/${testid}.dmp2  2>>$LOGFILE

  # !!15
  # dmlist $savfile header > $OUTDIR/${testid}.dmp1_std  2>>$LOGFILE
  # keyfilter $OUTDIR/${testid}.dmp1_std $OUTDIR/${testid}.dmp2_std \
  #            2>>$LOGFILE

  # compare
  # !!16
#   dmdiff $outfile $savfile tol=$SAVDIR/tolerance verb=0 > \
#         /dev/null 2>>$LOGFILE
#   if  test $? -ne 0 ; then
#     echo "ERROR: HEADER MISMATCH in $outfile" >> $LOGFILE
#     mismatch=0
#   fi

  ######################################################################
  # ascii files
  # !!17
  # diff $OUTDIR/${testid}.txt $OUTDIR/${testid}.txt_std > \
  #       /dev/null 2>>$LOGFILE
  # if  test $? -ne 0 ; then
  #   echo "ERROR: TEXT MISMATCH in $OUTDIR/${testid}.txt" >> $LOGFILE
  #   mismatch=0
  # fi

  ####################################################################
  # Did we get an error?
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

done
# end per-test loop
######################################################################


######################################################################
# report results

# blank line
echo ""

if test $script_succeeded -eq 0; then
    echo "${toolname} : PASS" | tee -a $LOGFILE
else
    echo "${toolname} : FAIL" | tee -a $LOGFILE
fi

echo "log file in ${LOGFILE}"


exit $script_succeeded
