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
# dmellipse.t
# test script for dmellipse


# !!2
# syntax:
# dmellipse.t [<testid> ... ]
 



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
toolname="dmellipse"

# set up list of tests
# !!4
alltests="test_one test_grid test_fix test_ang1 test_ang2 test_ang3 test_ang4 test_ang5 test_desc test_table_normal test_table_nullvals"

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

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do
    
  # delete old outputs
  rm -f $OUTDIR/${testid}*

  # Set up file names
  outfile=$OUTDIR/${testid}.reg
  savfile=$SAVDIR/${testid}.reg

  echo "running $testid" >> $LOGFILE

  ####################################################################
  # run the tool
  case ${testid} in
    # !!6
    test_one ) test1_string="dmellipse $INDIR/acisf00635_000N001_0001b_psf3.fits $outfile 0.9"
            ;;

    test_grid ) test1_string="dmellipse $INDIR/acisf00635_000N001_0001b_psf3.fits $outfile \"lgrid(0.5:0.91:0.05)\""
            ;;

    test_fix ) test1_string="dmellipse $INDIR/acisf00635_000N001_0001b_psf3.fits $outfile \"lgrid(0.5:0.91:0.05)\" fix_cen+ fix_ell+ fix_ang+"
            ;;

#
# This test checks the logic when the angle internally is negative
#

    test_ang1 ) test1_string="dmellipse $INDIR/hrcf01912_000N001_0002b_psf3.fits $outfile 0.95 norm- "
            ;;

    test_ang2 ) test1_string="dmellipse $INDIR/hrcf01912_000N001_0020b_psf3.fits $outfile 0.95 norm- "
            ;;

    test_ang3 ) test1_string="dmellipse $INDIR/hrcf01912_000N001_0024b_psf3.fits $outfile 0.95 norm- "
            ;;

    test_ang4 ) test1_string="dmellipse $INDIR/hrcf01912_000N001_0045b_psf3.fits $outfile 0.95 norm- "
            ;;

    test_ang5 ) test1_string="dmellipse $INDIR/hrcf01912_000N001_0151b_psf3.fits $outfile 0.95 norm- "
            ;;

    test_desc ) test1_string="dmellipse $INDIR/dmellipse_desc.fits $outfile 0.5"
            ;;

#
# KJG New for table stuff -- change Tes/INDIR to INDIR 
#
    test_table_normal ) test1_string="dmellipse Test/INDIR/'acisf00635_000N001_0001b_psf3_nonzero.tab[cols x,y,z=value]' $outfile 0.9 norm+ clob+"
            ;;
    
    test_table_nullvals ) test1_string="dmellipse Test/INDIR/i_mmm.tab.fits $outfile 0.9 norm+ clob+"
            ;;
    

  esac
  echo ""
  echo $test1_string | tee -a  $LOGFILE 
  echo ""
  eval $test1_string  | tee -a  $LOGFILE  2>&1
 

  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1

  # if different tests need different kinds of comparisons, use a 
  #  case ${testid} in...  here



dmdiff ${outfile} ${savfile}  tol=$SAVDIR/tolerance > /dev/null 2>>$LOGFILE
if  test $? -ne 0 ; then
  echo "ERROR: MISMATCH in $outfile" >> $LOGFILE
  mismatch=0
fi



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
