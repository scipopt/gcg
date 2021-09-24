#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#* This program is free software; you can redistribute it and/or             *
#* modify it under the terms of the GNU Lesser General Public License        *
#* as published by the Free Software Foundation; either version 3            *
#* of the License, or (at your option) any later version.                    *
#*                                                                           *
#* This program is distributed in the hope that it will be useful,           *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#* GNU Lesser General Public License for more details.                       *
#*                                                                           *
#* You should have received a copy of the GNU Lesser General Public License  *
#* along with this program; if not, write to the Free Software               *
#* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# @author Tim Donkiewicz
#
# This script generates a PDF with visualizations that compare two test runs
# that have been carried out using GCG's make test. It is called by the 
# Makefile, if make is executed with target "visu" and the data directory 
# contains a single run or the user chooses a single one.
# Alternatively, it is called via the check.sh script if make test
# is called with VISU=true.
#
# Common input arguments
# 1. file containing settings for this script (e.g. which visu to generate)
# 2-8. Test run statistics
#
# Additional input arguments from Makefile (if called via make visu)
# 9. Data directory
#
# Additional input arguments from check.sh (if called via make test VISU=true)
# 9.  path .out-file (dir and name), relative from inside check/
# 10. path to .res-file (dir and name), relative from inside check/
# 11. path to vbc dir, relative from inside check/
# 12-15. names and limits, given automatically by make test, calculated from runtime data in make visu

DEBUG="false" # enable extended logging for this script.

##############################################
####### Set variables, folders, names  #######
##############################################
# Current folder: check/

# Script settings file may include the variables below
SCRIPTSETTINGSFILE=$1
# Variables given by make test and make visu
BINID=$2
VERSION=$3
MODE=$4
LPS=$5
THREADS=$6
FEASTOL=$7
LAST_STATISTICS=$8
# variables only given by make test (will be calculated later if empty)
OUTFILE=$9
RESFILE=${10}
VBCFILES=${11}
TSTNAME=${12}
SETNAME=${13}
TIMELIMIT=${14}
MEMLIMIT=${15}

TIMESTAMP=$(date '+%d-%m-%Y_%H-%M-%S')
SCRIPTMODE=""

function echoGrey(){
  echo "$(tput dim)$1$(tput sgr0)"
}

function extractScriptSettings(){
  # This function will extract settings from a script settings file.
  # Though not recommended, you can also overwrite variables given by make test, e.g. to define your
  # own outfile/resfile/vbcfiles combination.
  #
  # Possible settings: (also set defaults before reading in settings file)
  ## Plot Toggles
  TREE=true       # generate tree plots?
  TIME=true       # generate time distribution plots?
  DETECTION=true  # generate detection plots?
  BOUNDS=true     # generate bounds plots?
  PRICING=true    # generate pricing plots?
  ## Plot Settings
  # for each plotter, give a flag that can contain arguments
  TREEARGS=
  TIMEARGS=
  DETECTIONARGS=
  BOUNDSARGS=
  PRICINGARGS=
  ## Proprietary plot settings (detection)
  CUSTOMCLASSIFIER=nonzeros     # for detection visus
  CUSTOMDETECTOR=SetPartMaster  # for detection visus (currently not configurable, fixed to SetPartMaster)
  ## General Settings
  DRAFT=false
  SET=$(realpath $SCRIPTSETTINGSFILE) # settings file has to lie inside check/

  # Check if file was given at all and if it exists
  if [ "$SCRIPTSETTINGSFILE" == "none" ]; then
    echo " No script settings file given. Using default script settings."
  elif test -f "$SET"; then
    echo "Script settings file $SCRIPTSETTINGSFILE found. Applying the following script settings:"
    # list changes
    awk -F\= '{gsub(/"/,"",$2);print "  " $1 " = " $2}' $SET
    # (re)define variables
    source $SET 
  else 
    echo " Could not find file $SCRIPTSETTINGSFILE. Using default script settings."
  fi

}

function setScriptMode(){
  # out file given
  # inside first if cond., because if given, we could be in make test VISU=true mode
  if [[ ! -z $OUTFILE && ($OUTFILE != "none") && (! -d $OUTFILE) ]]; then
    SCRIPTMODE="outfile"
  elif [[ ! -z $DATADIR ]]; then
    SCRIPTMODE="datadir"
  else
    if [[ (! -z $OUTFILE) && (! -z $DATADIR) ]]; then echo "No data given (or both, one run (out, res, vbc) AND data directory defined - dropping both)."; fi
    SCRIPTMODE="nodata"
  fi
}

function prepareVisualizationGeneration(){
  # function to get names of settings and testset from runtime files
  function getTstSetNames(){
    # based on runtime files, set test name and settings name if not yet set
    TSTNAME=$(cat $OUTFILE | sed -e 's/^M//g' -e 's/\r//g' | grep "saved parameter file" -m1 | cut -d "." -f 2 | sed "s/\_/\\\_/g")
    if [[ $TSTNAME = "" ]]; then TSTNAME=$(cat $OUTFILE | grep "saved parameter file" | cut -d '.' -f2 | head -n 1 | sed "s/\_/\\\_/g"); fi
    SETNAME=$(cat $OUTFILE | sed 's/^M//' | grep "loaded parameter file" -m1 | sed -e "s/>//g" | rev | cut -d '/' -f 1 | rev | sed -e "s/\.set//g" | sed "s/\_/\\\_/g")
    # if no settings file was read in (successfully! otherwise we could just use the filename, but it also prints settings if they don't exist), use default
    if [[ -z $SETNAME ]]; then SETNAME="default"; fi
  }

  # function to check a given data directory
  function checkDatadir(){
    # check if path is given and valid
    if [[ $DATADIR = "none" || $SCRIPTMODE = "nodata" ]]; then
        echo " No data directory given. The log directory should contain .res, .out and .vbc files for one or multiple runs."
        SCRIPTMODE=""
        return 1
    elif [[ ! -d $DATADIR ]]; then
        echo " Data directory could not be found."
        return 1
    fi

    DATADIR=$(realpath $DATADIR)
    # get number of files in data dir
    NOUTFILES=$(ls 2>/dev/null -Ubad1 -- $DATADIR/*.out | wc -l)
    NRESFILES=$(ls 2>/dev/null -Ubad1 -- $DATADIR/*.res | wc -l)
    NVBCFILES=$(ls 2>/dev/null -Ubad1 -- $DATADIR/*.vbc | wc -l)

    # check number of files in data dir
    if [[ $NOUTFILES = 0 ]]; then
        printf " Log directory does not contain any .res files. "
        read -p "Continue anyway? [Y/n] " yn
        if [[ $yn = "Y" ]]; then return 0; else return 1; fi
    elif test $NRESFILES = 0; then
        printf " Log directory does not contain any .out files. "
        read -p "Continue anyway? [Y/n] " yn
        if [[ $yn = "Y" ]]; then return 0; else return 1; fi
    else
        return 0
    fi
  }

  # function to check a given set of out file, res file and vbc files directory
  function checkOutRes(){
    # check if out and res files are given and valid
    if [ ! -f $OUTFILE ];   then read -e -p "  .out-file could not be found. Enter path to .out-file: " OUTFILE;
                            else OUTFILE=$(realpath $OUTFILE); fi
    if [ ! -f $RESFILE ];   then read -e -p "  .res-file could not be found. Enter path to .res-file: " RESFILE
                            else RESFILE=$(realpath $RESFILE); fi
    if [ ! -d $VBCFILES ];  then read -e -p "  .vbc-file could not be found. Enter path to vbc files directory: " VBCFILES
                            else VBCFILES=$(realpath $VBCFILES); fi
    getTstSetNames;
  }

  # function to check the vbc dirrectory, also trying commonly used ones
  # note: we require the standard vbc file name format (with settings name)
  # note: vbcdir is datadir by default (e.g. check/results/, completed by the function)
  function checkVBCdir(){
    VBCTRY=$VBCFILES
    # if no vbc files in vbcdir, try vbcdir/vbc/
    if test $(ls $VBCTRY/*.$SETNAME.vbc | wc -l) = 0; then
      VBCTRY=${VBCFILES}/vbc/; else VBCFILES=$VBCTRY; return 0
    fi
    # if no vbc files in vbcdir/vbc/, try vbcdir/settings/
    if test $(ls $VBCTRY/*.$SETNAME.vbc | wc -l) = 0; then
      VBCTRY=${VBCFILES}/$TSTNAME/; else VBCFILES=$VBCTRY; return 0
    fi
    # if no vbc files in vbcdir/settings/, try vbcdir/vbc/settings/ (default for GCG's testing)
    if test $(ls $VBCTRY/*.$SETNAME.vbc | wc -l) = 0; then
      VBCTRY=${VBCFILES}/vbc/$TSTNAME/; else VBCFILES=$VBCTRY; return 0
    fi
    # no more possibilities to test
    if test $(ls $VBCTRY/*.$SETNAME.vbc | wc -l) = 0; then
      VBCFILES=""; return 1; else VBCFILES=$VBCTRY; return 0
    fi
  }

  # this function checks if in the data directory, there is only one run
  # prompting to choose one if there is more than one run
  function checkNRuns(){
    VBCFILES_def=$VBCFILES # save user input, if applicable
    outfiles=($(ls $(realpath $DATADIR)/*.out))
    if [[ ($NOUTFILES > 1) && ($NRESFILES > 1) ]]; then
      i=1
      echo "  Data directory contains more than one run:"
      for file in ${outfiles[@]}; do
        VBCFILES=$VBCFILES_def # reinitialize as user input, if applicable
        # initialize runtime data variables
        OUTFILE=$file
        RESFILE=${OUTFILE%.out}.res # resfile must be named like outfile
        getTstSetNames;
        if [ -z $VBCFILES_def ]; then VBCFILES=$DATADIR; checkVBCdir > /dev/null 2>&1; else VBCFILES=VBCFILES_def; checkVBCdir; fi
        # print them to check
        echo "  [$i]" $(basename $OUTFILE)
        if [[ -f $RESFILE ]]; then
          echoGrey "      $(basename $(realpath $RESFILE))"; else
          echoGrey "      (no resfile)"
        fi
        if [[ -d $VBCFILES ]]; then
          echoGrey "      $(realpath $VBCFILES)"; else
          echoGrey "      (no valid vbc dir)"
        fi
        ((i=i+1))
      done
      read -p "  Choose a number for testset report or type 'c' to generate comparison report. " n
      if [[ $n = 'c' ]]; then
        ./writeComparisonReport.sh $SCRIPTSETTINGSFILE $LAST_STATISTICS $DATADIR
        exit 0
      else 
        OUTFILE=${outfiles[(($n-1))]}
        RESFILE=${OUTFILE%.out}.res # resfile must be named like outfile
        if [ -z $VBCFILES ]; then VBCFILES=$DATADIR; fi
      fi
    else
        # if there are multiple res files but just one out file
        OUTFILE=${outfiles[0]}
        RESFILE=${OUTFILE%.out}.res # resfile must be named like outfile
        if [ -z $VBCFILES ]; then VBCFILES=$DATADIR; fi
    fi
  }

  # check if the given runtime data can be found
  # default: if none given, ask for data directory
  setScriptMode
  if [[ $SCRIPTMODE = "datadir" || $SCRIPTMODE = "nodata" ]]; then
    nchecks=0
    checkDatadir
    while [ $? -eq 1 ]; do
      # check if data directory was entered incorrectly five times, if so, terminate
      if [ $nchecks = 5 ]; then
        echo "  No valid directory was given five times in a row. Terminating."
        exit 1
      else
        ((nchecks=nchecks+1))
      fi
      # else try to read it again
      read -e -p "  Enter path to log files directory: " DATADIR;
      checkDatadir
    done
    # check the number of runs. if more then two are found, suggest comparison report
    checkNRuns
  elif [[ $SCRIPTMODE = "outfile" ]]; then
    checkOutRes
    checkVBCdir
  fi

  # get testset and settings names to then set report directory 
  getTstSetNames
  if [ -z $REPORTDIR ]; then
    mkdir -p reports/ > /dev/null 2>&1
    REPORTDIR=$(realpath "reports/testsetreport_${TSTNAME}_${SETNAME}_${TIMESTAMP}") # folder to put report into (default: inside check/reports)
    mkdir -p $REPORTDIR > /dev/null 2>&1
    printf " Report output folder: ${REPORTDIR} \n"
  fi

  # check existence of vbc files (only for test set report)
  if test $DEBUG = "true"; then
    checkVBCdir
  else
    checkVBCdir > /dev/null 2>&1
  fi
  while [ $? -eq 1 ]; do
    echo "  Invalid vbc directory or no .vbc files in it."
    read -e -p "  Enter path to vbc files directory (or press return to skip - last dir will be used): " VBCFILES;
    if [[ $VBCFILES != "" ]]; then 
      if test $DEBUG = "true"; then
        checkVBCdir
      else
        checkVBCdir > /dev/null 2>&1
      fi
    else
      NOVBC=true
      break
    fi
  done
  
  # define output folders and file
  PLOTDIR=${REPORTDIR}/plots
  LOGDIR=${REPORTDIR}/logs
  VBCDIR=${LOGDIR}/vbc
  REPORTFILE=${TSTNAME}.${SETNAME}.testsetreport.tex

  # create output folders
  mkdir -p ${PLOTDIR} > /dev/null 2>&1
  mkdir -p ${LOGDIR} > /dev/null 2>&1
  if [[ -z $NOVBC ]]; then mkdir -p ${VBCDIR} > /dev/null 2>&1; fi

  # make paths absolute
  OUTFILE=$(realpath $OUTFILE)
  RESFILE=$(realpath $RESFILE)
  if [[ -z $NOVBC ]]; then VBCFILES=$(realpath $VBCFILES); fi

  # move runtime data files to log directory inside report folder
  cp ${OUTFILE} ${LOGDIR} -r
  cp ${RESFILE} ${LOGDIR} -r
  if [[ -z $NOVBC ]]; then cp ${VBCFILES}/*.vbc ${VBCDIR} -r; fi

  # after moving, we do not need to access check/results anymore
  OUTFILE=$(ls ${LOGDIR}/*.out)
  RESFILE=$(ls ${LOGDIR}/*.res)

  # Variables to generate for the front page of the report
  # get general limits from outfile
  TIMELIMIT=$(cat $OUTFILE | grep "set limits time" -m1 | cut -d " " -f 5)
  MEMLIMIT=$(cat $OUTFILE | grep "set limits memory" -m1 | cut -d " " -f 5)
  NODELIMIT=$(cat $OUTFILE | grep "set limits nodes" -m1 | cut -d " " -f 5)
  # get overall statistics from resfile
  INSTANCES=$(sed -e '/^--/d' -e '/^Name/d' -e "/@/d" ${RESFILE} | cut -d' ' -f1 | sed '/^$/d')
  NINSTANCES=$(echo $(tail -n2 $RESFILE) | cut -d " " -f 1)
  NINSTANCES=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 1)
  NINSTANCES_PASS=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 2)
  NINSTANCES_LIMIT=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 3)
  NINSTANCES_FAIL=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 4)
  TOTSOLVE=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 7)
  AVGSOLVE=$(echo $(tail -n5 $RESFILE | head -n1) | cut -d " " -f 8)
}

function generateVisualizations(){
##############################################
# Plot everything that's not comparing runs: #
##############################################
# visu scripts must always be executed from within the stats folder.
cd ../stats
# Current folder: stats/
  if test $DRAFT = "true"; then
    echo " [draft mode: time bar chart and bounds plot only]"
      echo -ne '|░░░░░               |  (25%)  Time Plot     \r'
      python3 general/time.py ${OUTFILE} --outdir $PLOTDIR/timedist --bar $TIMEARGS > /dev/null 2>&1
      echo -ne '|░░░░░░░░░░░░░░░     |  (75%)  Bounds Plot    \r'
      python3 bounds/plotter_bounds.py ${OUTFILE} --outdir ${PLOTDIR}/bounds $BOUNDSARGS > /dev/null 2>&1
      echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done generating visualizations     \r'
      echo -ne '\n'
  elif test $DEBUG = "true"; then
    if test $LAST_STATISTICS = "false"; then
      echo " [debug mode: no compilation with STATISTICS=true, skipping bounds and pricing visualizations]"
    else
      echo " [debug mode: all]"
    fi
    if test $TREE = "true"; then echo 'Tree Plot'
    python3 tree/plotter_tree.py ${VBCDIR} --outdir ${PLOTDIR}/tree $TREEARGS ; fi
    if test $TIME = "true"; then echo 'Time Plot'
    python3 general/time.py ${OUTFILE} --outdir $PLOTDIR/timedist -A $TIMEARGS;
    python3 general/time.py ${OUTFILE} --outdir $PLOTDIR/timedist --single $TIMEARGS; fi
    if test $DETECTION = "true"; then echo 'Detection Plot'
    python3 detection/plotter_detection.py ${OUTFILE} --outdir ${PLOTDIR}/detection $DETECTIONARGS; fi
    # The following scripts require GCG to be compiled with STATISTICS=true
    if test $LAST_STATISTICS = "true"; then
      if test $BOUNDS = "true"; then echo 'Bounds Plot'
      python3 bounds/plotter_bounds.py ${OUTFILE} --outdir ${PLOTDIR}/bounds $BOUNDSARGS;
      python3 bounds/plotter_bounds.py ${OUTFILE} --outdir ${PLOTDIR}/bounds --xaxis=iter $BOUNDSARGS; fi
      if test $PRICING = "true"; then echo 'Pricing Plot'
      python3 pricing/plotter_pricing.py ${OUTFILE} --no-vartime --vbcdir=${VBCDIR} --outdir ${PLOTDIR}/pricing $PRICINGARGS; fi
    fi
  else
    if test $LAST_STATISTICS = "false"; then
      echo " [no compilation with STATISTICS=true, skipping bounds and pricing visualizations]"
    else
      echo " [all]"
    fi
    if test $TREE = "true"; then echo -ne '|░░                  |  (10%)  Tree Plot     \r'
    python3 tree/plotter_tree.py ${VBCDIR} --outdir ${PLOTDIR}/tree $TREEARGS > /dev/null 2>&1; fi
    if test $TIME = "true"; then echo -ne '|░░░░                |  (20%)  Time Plot     \r'
    python3 general/time.py ${OUTFILE} --outdir $PLOTDIR/timedist -A $TIMEARGS > /dev/null 2>&1; fi
    if test $TIME = "true"; then echo -ne '|░░░░░░              |  (30%)  Time Plot     \r'
    python3 general/time.py ${OUTFILE} --outdir $PLOTDIR/timedist --single $TIMEARGS > /dev/null 2>&1; fi
    if test $DETECTION = "true"; then echo -ne '|░░░░░░░░            |  (40%)  Detection Plot   \r'
    python3 detection/plotter_detection.py ${OUTFILE} --outdir ${PLOTDIR}/detection -c $CUSTOMCLASSIFIER $DETECTIONARGS> /dev/null 2>&1; fi
    # The following scripts require GCG to be compiled with STATISTICS=true
    if test $LAST_STATISTICS = "true"; then
      if test $BOUNDS = "true"; then echo -ne '|░░░░░░░░░░          |  (50%)  Bounds Plot    \r'
      python3 bounds/plotter_bounds.py ${OUTFILE} --outdir ${PLOTDIR}/bounds $BOUNDSARGS > /dev/null 2>&1; fi
      if test $BOUNDS = "true"; then echo -ne '|░░░░░░░░░░░░        |  (60%)  Bounds Plot    \r'
      python3 bounds/plotter_bounds.py ${OUTFILE} --outdir ${PLOTDIR}/bounds --xaxis=iter $BOUNDSARGS > /dev/null 2>&1; fi
      if test $PRICING = "true"; then echo -ne '|░░░░░░░░░░░░░░      |  (70%)  Pricing Plot   \r'
      python3 pricing/plotter_pricing.py ${OUTFILE} --no-vartime --vbcdir=${VBCDIR} \
      --outdir ${PLOTDIR}/pricing $PRICINGARGS > /dev/null 2>&1; fi
    fi
    echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done generating visualizations     \r'
    echo -ne '\n'
  fi
}

function printTeX_makefile(){
cat > ${TSTNAME}.${SETNAME}.testsetreport.make << EndOfMessage

# latexmk automatically manages the .tex files
${TSTNAME}.${SETNAME}.testsetreport.pdf: ${REPORTFILE}
	@echo ------------
	@echo 
	@echo Compiling tex code. This may take a while.
	@echo 
	@echo ------------
	@latexmk -f -pdf -pdflatex="pdflatex -interaction=batchmode -shell-escape" -use-make ${REPORTFILE}
	@make -f ${TSTNAME}.${SETNAME}.testsetreport.make clean

clean:
	@latexmk -c
	@rm -f report_*figure*.*
	@rm -f *.auxlock
	@rm -f *figure*.md5
	@rm -f *figure*.log
	@rm -f *figure*.dpth

cleanall:
	@latexmk -C
	@make -f ${TSTNAME}.${SETNAME}.testsetreport.make clean
EndOfMessage
}

function printTeX_readme(){
cat > ${TSTNAME}.${SETNAME}.testsetreport.README << EndOfMessage
README: How to create a PDF file from the .tex file(s) using the ${TSTNAME}.make file.
Note: The package pdflatex is required.

Use the command
    'make -f ${TSTNAME}.${SETNAME}.testsetreport.make'
to compile.
Depending on the size of your problem that may take some time.
Please do not delete any new files that might be generated during the compile process.
All access files will be deleted automatically once the compilation is complete.

Clean options:
    'make -f ${TSTNAME}.${SETNAME}.testsetreport.make clean' clears all present intermediate files (if any exist)
    'make -f ${TSTNAME}.${SETNAME}.testsetreport.make cleanall' clears all generated files INCLUDING .pdf
EndOfMessage
}

function printTeX_preamble(){
cat > ${REPORTFILE} << EndOfMessage
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                           *
% *                  This file is part of the program                         *
% *          GCG --- Generic Column Generation                                *
% *                  a Dantzig-Wolfe decomposition based extension            *
% *                  of the branch-cut-and-price framework                    *
% *         SCIP --- Solving Constraint Integer Programs                      *
% *                                                                           *
% * Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       *
% *                         Zuse Institute Berlin (ZIB)                       *
% *                                                                           *
% * This program is free software; you can redistribute it and/or             *
% * modify it under the terms of the GNU Lesser General Public License        *
% * as published by the Free Software Foundation; either version 3            *
% * of the License, or (at your option) any later version.                    *
% *                                                                           *
% * This program is distributed in the hope that it will be useful,           *
% * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
% * GNU Lesser General Public License for more details.                       *
% *                                                                           *
% * You should have received a copy of the GNU Lesser General Public License  *
% * along with this program; if not, write to the Free Software               *
% * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*
% *                                                                           *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%
% @author Tim Donkiewicz

\documentclass[a4paper,10pt]{article}

% packages
\usepackage[utf8]{inputenc}
\usepackage[hidelinks]{hyperref}
\usepackage{pdfpages}
\usepackage{fancybox}
\usepackage{float}
\setcounter{secnumdepth}{4}
\usepackage{graphicx}
\usepackage{grffile}

\begin{document}

\begin{titlepage}
  \centering
  \thispagestyle{empty}
  {\Huge Test Set Report} \\\\ \today

\vspace{2cm}
\begin{tabular}{{lp{10cm}}}
    \multicolumn{2}{l}{\large\textbf{Test Run Characteristics}} \\\\
    \textbf{Test Set Name}: & \begin{minipage}{10cm}\begin{verbatim}${TSTNAME}\end{verbatim}\end{minipage} \\\\
    \textbf{Settings}:      & \begin{minipage}{10cm}\begin{verbatim}${SETNAME}\end{verbatim}\end{minipage} \\\\
    Mode:                   & \begin{minipage}{10cm}\begin{verbatim}${MODE}\end{verbatim}\end{minipage} \\\\
    GCG Version:            & \begin{minipage}{10cm}\begin{verbatim}${VERSION}\end{verbatim}\end{minipage} \\\\
    Feasibility Tolerance:  & ${FEASTOL} \\\\
    LP Solver:              & ${LPS} \\\\
    Time Limit:             & ${TIMELIMIT}s\\\\
    Memory Limit:           & ${MEMLIMIT}MB \\\\
    \vspace{0.5cm}
    Node Limit:             & ${NODELIMIT}\\\\

    \multicolumn{2}{l}{\large\textbf{Test Set Solving Overview}} \\\\
    Number of Instances:    & ${NINSTANCES} (${NINSTANCES_PASS} passed, ${NINSTANCES_FAIL} failed, ${NINSTANCES_LIMIT} ran into a limit)\\\\
    Total Solving Time:     & ${TOTSOLVE}s \\\\
    \vspace{1cm}
    Average Solving Time:   & ${AVGSOLVE}s \\\\

    Report Timestamp:       & $(echo ${TIMESTAMP} | sed "s/\_/ /g" )\\\\
    Test Run Binary ID:     & $(echo ${BINID} | sed "s/\_/\\\_/g" ) \\\\
    
\end{tabular}

\end{titlepage}
\newpage
\thispagestyle{empty}
\tableofcontents
\newpage
EndOfMessage
}

# define descriptions for each visualization (deactivate automatic line breaks in your IDE for improved readibility)
declare -A desc_inst
desc_inst["bounds-time"]="The top subplot shows the development of the primal and dual bounds in the RMP during the pricing in the root node as given by the table ``root bounds'' printed by GCG. Every change represents a pricing iteration and the resulting changes to the bounds. The bounds are complemented by a newly created gap plot, which will be explained in Section {sec:tgpp}. The other two subplots illustrate the point in time in the pricing at which the columns that are finally in the basis are generated."
desc_inst["bounds-iter"]="The same plot as the bounds time plot (see below), but with pricing iterations in the root node instead of the time spent there on the \$x\$-axis."
desc_inst["tree-bar"]="This plot shows the percentage of nodes in the Branch-and-Bound tree opened on each level against how many exist on this level (\$2^x\$)."
desc_inst["tree-plot"]="This plot shows the distribution of nodes in the Branch-and-Bound tree opened in absolute terms."
desc_inst["pricing-complete"]="This Plot shows how many variables were generated in a certain pricing round in which time for all nodes of the Branch and Bound tree. The node numbers are shown above the plot and the rounds are in the line below that. Each bar represents the iteration of one pricing problem. Note that the numbers of the pricing problems can have gaps in between, since they could have been aggregated prior to the pricing. Whether those variables are useful is shown by all bars that are below zero, as they mean that the variables of that pricing iteration are in the optimal solution of the Root LP (Root LP Sol) or IP (Incumbent). Finally, the dots show how many columns are taken from the column pool."
desc_inst["pricing-time"]="The Pricing Time Statistics include four pie charts. The first one shows how much of the runtime was needed in the reduced cost pricing, the master LP and during the initial Farkas. The upper center one shows the relative (and, inside the slices, absolute) time needed by each pricing problem that took at least \$\frac{11}{360}\$ of the total pricing time (\$11^{\circ}\$ of the pie, the last degree where the absolute numbers inside the slices are still readable). Note that if no absolute numbers are needed, but only the highest possible amount of slices (pricing problems) should be shown, the \texttt{-{}-short-times} argument can be set. The pie chart to the upper right shows how many columns were generated by each pricing problem and the ratio between the upper right and the upper center, i.e. the variables per second, is shown in the lower left, illustrating which pricing problem yielded the most variables for the RMP. Finally, in the course of this thesis, an additional subplot that illustrates the seconds needed by each pricing problem to generate a variable was added."
desc_inst["pricing-summary"]="The summary plot aims to illustrate the same thing as the ``complete plot''. The end of the root node, which is treated in deeper detail in the Bounds Plot, is marked by a red line. The plot consists of two different \$y\$-axes, one representing the time (in seconds) needed for the pricing and the other the fraction of pricing problems that generated variables. This leads to the ability to identify pricing rounds that ran for a long time and see when and how many pricing problems were successful."
desc_inst["pricing-bubble"]="In this visualization, one can see all pricing problems listed vertically along the \$y\$-axis. Then, in the left subfigure, they are shown against the pricing rounds on the \$x\$-axis. Every time the pricer yielded at least one variable resulting from a pricing problem, a dot is printed in the round where it was generated. This results in the ability to not only see the sensibility of each pricing problem, but also in which rounds what pricing problem performed best. The subplot on the right-hand side shows how many percent of the variables were generated by which problem."
desc_inst["pricing-gap"]="The gap plot compares the solving time of a pricing problem (of which usually more than one are solved during a pricing \textit{round}) with the size of the gap in the root node at that point of time. The gap here is the ratio between the maximum gap (the ratio of the worst, but finite primal bound and the worst, but finite dual bound) and the gap in the root node at the point of time given by the \$x\$-axis. Note that it is also possible to show the time of one pricing \textit{iteration} instead of one pricing \textit{problem}, which sums up all points seen here that happened in the same iteration. Note that the data used for the plots in general is only taken once per \$0.01\$ seconds, which leads to a ``discrete'' distribution."
desc_inst["pricing-depth"]="This figure illustrates how the gap develops along the depth of the branching tree. Each dot represents the gap as given by the primal and dual bounds in this specific node as given by the GCG ``root bounds'' table (just like in the bounds plotter). This node is located on the tree depth that can be read on the \$x\$-axis, such that for each \$x\$-coordinate, at most \$2^x\$ points can exist. Furthermore, a plot of the mean is given."
desc_inst["pricing-nodeID"]="The Node ID plot is similar to the Depth Plot. Instead of the depth in the branch-and-bound tree, we now have the node ID. This leads to the fact that one can see behavior that is not dependent of the depth, but of the time progression during the branching."

declare -A desc_agg
desc_agg["detection-times"]="This plot visualizes what fraction of instances against the time used for the whole detection process, including classification and score computation, on a logarithmic scale."
desc_agg["detection-decomps"]="The above plot shows the number of decompositions that were found for what fraction of instances in the given test set on a logarithmic scale. Only those decompositions are shown that have a score that is strictly greater than \$0\$."
desc_agg["detection-quality"]="This plot illustrates the detection goodness of the test set (according to the ``max white score''). If for all instances in the test set, the best (``whitest'') decomposition was completely white, this plot would show a line with \$y=1\$ for all \$x\$."
desc_agg["detection-quality_custom"]="Similarly to the previous plot, this one shows the detection goodness of the whole test set, but for a specific detector, the ``$CUSTOMDETECTOR'' detector. For this specific detector, GCG outputs the scores separately in the \texttt{detectionstatistics} test mode. With this plot, together with the previous one, one can compare the performance of the Set Partitioning Master detector with the overall performance."
desc_agg["detection-nBlocksOfBest"]="This plot shows how many blocks are used in the (according to the max white score) best decomposition, allowing to make statements about the sensibility of different numbers of blocks."
desc_agg["detection-classification_classes_custom"]="With this visualization, one can see how many classes a classifier determined for the variables. It can be chosen for which classifier to plot this, here ``$CUSTOMCLASSIFIER'' was chosen."
desc_agg["timedist-plot"]="A simple plot of normalized times used for each instance of the testset."
desc_agg["timedist-bar"]="A stacked bar chart of the normalized time distribution of all instances in the test set."
desc_agg["timedist-pie"]="An averaged pie chart for the whole testset."
desc_agg["timedist-grouped_bar"]="A normalized grouped bar chart of the testset, showing the distribution of the distribution of times."


function printTeX_instance_information(){
# start generating instance-related
# types and subtypes of visus:
# - bounds plot (2 per instance)
# - pricing plots (6 per instance)
# - tree plots (2 per instance)

# define which sections to create in TeX report
visus_inst=''
if [[ $TREE = "true" ]]; then visus_inst="$visus_inst tree"; fi
if [[ ($BOUNDS = "true") && ($LAST_STATISTICS = "true") ]]; then visus_inst="$visus_inst bounds"; fi
if [[ ($PRICING = "true") && ($LAST_STATISTICS = "true") ]]; then visus_inst="$visus_inst pricing"; fi

echo "\section{Instance Information}" >> ${REPORTFILE}
for i in ${INSTANCES} # instances
do
  echo "\subsection{Instance: $(echo ${i} | sed "s/\_/\\\_/g" )}" >> ${REPORTFILE}
 # echo "\addcontentsline{toc}{section}{Instance: $(echo ${i} | sed "s/\_/\\\_/g" )}" >> ${REPORTFILE}
  for t in $visus_inst # type of visu
  do
    for v in $(ls plots/$t/$i*) # subtype of visu
    do
      echo "\begin{figure}[H]"  >> ${REPORTFILE} 
      echo "\makebox[\textwidth][c]{\includegraphics[width=1.2\textwidth]{$v}}" >> ${REPORTFILE}
      echo "  \caption{${desc_inst["${t}-$(echo "$v" | rev | cut -d. -f2 | rev)"]} \\\\Visualization Path: \texttt{$(echo "$v" | sed "s/\_/\\\_/g")}}" >> ${REPORTFILE}
      echo "\end{figure}"       >> ${REPORTFILE}
      echo "\vspace*{\fill}\newpage\vspace*{\fill}" >> ${REPORTFILE}
    done
  done
done
}

function printTeX_aggregated_information(){
# start generating aggregated
# types and subtypes of visus:
# - time distribution (5 total)
# - detection statistics (6 total)

# define which sections to create in TeX report
visus_agg=''
if [[ $TIME = "true" ]]; then visus_agg="$visus_agg timedist"; fi
if [[ $DETECTION = "true" ]]; then visus_agg="$visus_agg detection"; fi

declare -A titles
titles["timedist"]="Time Distribution Plot"
titles["detection"]="Detection Visualizations"

echo "\newpage\section{Aggregated Information}" >> ${REPORTFILE}
for p in $visus_agg
do
  echo "\subsection{${titles["${p}"]}}" >> ${REPORTFILE}
  for v in $(ls plots/$p/* | sed "/singlebar/d") # subtype of visu
  do
    echo "\vspace*{\fill}"    >> ${REPORTFILE} 
    echo "\begin{figure}[H]"  >> ${REPORTFILE} 
    echo "\makebox[\textwidth][c]{\includegraphics[width=1.2\textwidth]{$v}}" >> ${REPORTFILE}
    echo "  \caption{${desc_agg["${p}-$(echo "$v" | rev | cut -d. -f2 | rev | sed -e "s/$CUSTOMDETECTOR/custom/g" -e "s/$CUSTOMCLASSIFIER/custom/g")"]}\\\\Visualization Path: \texttt{$(echo "$v" | sed "s/\_/\\\_/g")}}" >> ${REPORTFILE}
    echo "\end{figure}"       >> ${REPORTFILE}
    echo "\vspace*{\fill}\newpage" >> ${REPORTFILE}
  done
done
}

function printTeX_tail(){
  echo "\end{document}" >> ${REPORTFILE}
}


function generateReport(){
  cd ${REPORTDIR}
  echo -ne '|░░                  |  (10%)  Generating TeX Report Code... \r'
  printTeX_makefile; 
  printTeX_readme;
  printTeX_preamble;
  echo -ne '|░░░░                |  (20%)  Generating TeX Report Code... \r'
  if [ $NINSTANCES -gt 1 ]; then
    printTeX_aggregated_information > /dev/null 2>&1
  fi
  echo -ne '|░░░░░░░░            |  (40%)  Generating TeX Report Code... \r'
  # if folders are empty, there was no detection/... for that instance! (we do not output empty folders)
  printTeX_instance_information > /dev/null 2>&1
  printTeX_tail;
  echo -ne '|░░░░░░░░░░░░        |  (60%)  Compiling TeX Report Code...  \r'
  make -f ${TSTNAME}.${SETNAME}.testsetreport.make > /dev/null 2>&1
  echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done! Opening Report...       \r'
  echo -ne '\n'
  xdg-open ${TSTNAME}.${SETNAME}.testsetreport.pdf > /dev/null 2>&1
  if [ $? -eq 1 ]; then echo "Opening PDF.";
  else echo "Could not open PDF file."; fi
}


function main(){
  printf "\n$(tput bold)Starting Testset Report Generation...$(tput sgr0)\n"
  extractScriptSettings;
  prepareVisualizationGeneration;
  printf "\n$(tput bold)Generating visualizations.$(tput sgr0)"
  generateVisualizations;
  printf "\n$(tput bold)Generating PDF report file.\n$(tput sgr0)"
  generateReport;
}

main;