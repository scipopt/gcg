#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       *
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
# contains more than one run.
#
# Input arguments from Makefile:
# 1. file containing settings for this script (e.g. which visu to generate)
# 2. flag to indicate whether your current compile is with STATISTICS=true
#    (if not, it will make the script not generate bounds plots)
# 3. folder with .out and .res files to compare runs of

DEBUG="false" # enable extended logging for this script.

##############################################
####### Set variables, folders, names  #######
##############################################
# Current folder: check/

# Script settings file may include the variables below
SCRIPTSETTINGSFILE=$1
LAST_STATISTICS=$2
DATADIR=$3

TIMESTAMP=$(date '+%d-%m-%Y_%H-%M-%S')

function echoBold(){
  if [ "$DEBUG" == "false" ]; then printf "\n$(tput bold)$1$(tput sgr0)"
  else printf "\n$1"; fi
}

function extractScriptSettings(){
  # This function will extract settings from a script settings file.
  # Though not recommended, you can also overwrite variables given by make test, e.g. to define your
  # own outfile/resfile/vbcfiles combination.
  #
  # Possible settings: (also set defaults before reading in settings file)
  ## Plot Toggles
  TABLE=true      # generate comparison table?
  GENERAL=true    # generate general plots?
  PERFPROF=true   # generate performance profile?
  TIME=true       # generate time distribution plots?
  DETECTION=true  # generate detection plots?
  BOUNDS=true     # generate bounds plots?
  ## Plot Settings
  # for each plotter, give a flag that can contain arguments
  PERFPROFARGS=
  TIMEARGS=
  DETECTIONARGS=
  DETECTIONARGS=
  BOUNDSARGS=
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
    echo " Script settings file $SCRIPTSETTINGSFILE found. Applying the following script settings:"
    # list changes
    awk -F\= '{gsub(/"/,"",$2);print "  " $1 " = " $2}' $SET
    # (re)define variables
    source $SET
  else
    echo " Could not find file $SCRIPTSETTINGSFILE. Using default script settings."
  fi
}

function prepareVisualizationGeneration(){
  # get runtime files from user
  function checkDatadir(){
    if [ -d $DATADIR ]; then DATADIR=$(realpath $DATADIR);
    elif [ -d ../$DATADIR ]; then DATADIR=$(realpath ../$DATADIR); fi

    if [[ $DATADIR = "none" ]]; then
        echo " No log directory given. The log directory should contain .res and .out files for all runs."
        return 1
    elif [[ ! -d $DATADIR ]]; then
        echo " Log directory could not be found."
        return 1
    elif [[ $(ls 2>/dev/null -Ubad1 -- $(realpath $DATADIR)/*.res | wc -l) = 0 ]]; then
        printf " Log directory does not contain any .res files. "
        read -p "Continue anyway? [Y/n] " yn
        if [[ $yn = "Y" ]]; then return 0; else return 1; fi
    elif test $(ls 2>/dev/null -Ubad1 -- $(realpath $DATADIR)/*.out | wc -l) = 0; then
        printf " Log directory does not contain any .out files. "
        read -p "Continue anyway? [Y/n] " yn
        if [[ $yn = "Y" ]]; then return 0; else return 1; fi
    else
        return 0
    fi
  }
  # check if the given log directory is valid, else ask for one
  checkDatadir
  while [ $? -eq 1 ]; do
    read -e -p "  Enter path to log files directory: " DATADIR;
    checkDatadir
  done

  # get testset and settings names to set report directory
  if [ -z $REPORTDIR ]; then
    mkdir -p reports/ > /dev/null 2>&1
    REPORTDIR=$(realpath "reports/comparisonreport_${TIMESTAMP}") # folder to put report into (default: inside check/reports)
    mkdir -p $REPORTDIR > /dev/null 2>&1
    printf " Report output folder: ${REPORTDIR} \n"
  fi

  # define output folders and file
  DATADIR=$(realpath $DATADIR)
  PLOTDIR=$(realpath ${REPORTDIR}/plots)
  LOGDIR=$(realpath ${REPORTDIR}/logs)
  PKLDIR=$(realpath ${REPORTDIR}/pickles)
  REPORTFILE=comparisonreport.tex

  # create output folders
  mkdir -p ${PLOTDIR} > /dev/null 2>&1
  mkdir -p ${LOGDIR} > /dev/null 2>&1
  mkdir -p ${PKLDIR} > /dev/null 2>&1

  # move runtime data files to log directory inside report folder
  cp ${DATADIR}/*.res ${LOGDIR} -r  > /dev/null 2>&1
  cp ${DATADIR}/*.out ${LOGDIR} -r  > /dev/null 2>&1

  # after moving, we do not need to access check/results anymore
  OUTFILES=$(ls ${LOGDIR}/*.out)
  RESFILES=$(ls ${LOGDIR}/*.res)

  # Variables to generate for the front page of the report
  NRUNS=$(ls 2>/dev/null -Ubad1 -- $(realpath $DATADIR)/*.res | wc -l)
}


function generateVisualizations(){
    ##############################################
    # Plot everything that's not comparing runs: #
    ##############################################
    chmod +x parseres.py
    chmod +x plotcomparedres.py

    if test $DRAFT = "true"; then
        echo " [draft mode: performance profile only]"
        echo -ne '|░░░░░░░░░░          |  (50%)  Performance Profile \r'
        python3 ../stats/general/performance_profile.py $(ls $RESFILES) --outdir "$(realpath ${PLOTDIR})/" > /dev/null 2>&1
        echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done generating visualizations     \r'
        echo -ne '\n'
    elif test $DEBUG = "true"; then
        if test $LAST_STATISTICS = "false"; then
          echo " [debug mode: no compilation with STATISTICS=true, skipping bounds visualizations]"
        else
          echo " [debug mode: all]"
        fi
        if test $GENERAL = "true"; then echo 'General Plots'
            for i in $RESFILES; do ./parseres.py $i ${PKLDIR}; done
            python3 plotcomparedres.py ${PKLDIR} ${PLOTDIR}/general/;
        fi
        # visu scripts should be executed from within the stats folder.
        cd ../stats
        # generate comparison table (first remove all old tables, since this scripts output dir is static)
        if test $TABLE = "true"; then echo 'Comparison Table'
        rm -f plots/*.tex
        ./general/comparison_table.sh $RESFILES  > /dev/null 2>&1 && mv plots/*.tex $PLOTDIR; fi
        if test $PERFPROF = "true"; then echo 'Performance Profile'
        python3 general/performance_profile.py $RESFILES --outdir ${PLOTDIR}/performance_profile/ $PERFPROFARGS; fi
        if test $TIME = "true"; then echo 'Time Plot: Parsing'
        for j in $OUTFILES; do python3 general/time.py $j --outdir $PKLDIR --save $TIMEARGS; done; fi
        if test $TIME = "true"; then echo 'Time Plot: Plotting'
        python3 general/time.py $PKLDIR/check.*.pkl --outdir $PLOTDIR/timedist/ --compare $TIMEARGS; fi
        if test $TIME = "true"; then echo 'Time Plot: Plotting'
        python3 general/time.py $PKLDIR/check.*.pkl --outdir $PLOTDIR/timedist/ --compare $TIMEARGS; fi
        if test $DETECTION = "true"; then echo 'Detection Plot'
        python3 detection/plotter_detection.py $OUTFILES --outdir $PLOTDIR/detection/ $DETECTIONARGS; fi
        # The following scripts require GCG to be compiled with STATISTICS=true
        if test $LAST_STATISTICS = "true"; then
            if test $BOUNDS = "true"; then echo 'Bounds Plot'
            python3 bounds/plotter_bounds.py $OUTFILES --outdir ${PLOTDIR}/bounds $BOUNDSARGS; fi
        fi
        cd ../check
        echo 'Done generating visualizations'
    else
        if test $LAST_STATISTICS = "false"; then
          echo " [no compilation with STATISTICS=true, skipping bounds visualizations]"
        else
          echo " [all]"
        fi
        if test $GENERAL = "true"; then echo -ne '|░░                  |  (10%)  General Plots        \r'
            for i in $RESFILES; do ./parseres.py $i ${PKLDIR} > /dev/null 2>&1; done
            python3 plotcomparedres.py ${PKLDIR} ${PLOTDIR}/general/  > /dev/null 2>&1;
        fi
        # visu scripts should be executed from within the stats folder.
        cd ../stats
        # generate comparison table (first remove all old tables, since this scripts output dir is static)
        if test $TABLE = "true"; then echo -ne '|░░░░░░              |  (20%)  Comparison Table    \r'
        rm -f plots/*.tex
        ./general/comparison_table.sh $RESFILES  > /dev/null 2>&1 && mv plots/*.tex $PLOTDIR; fi
        if test $PERFPROF = "true"; then echo -ne '|░░░░░░              |  (30%)  Performance Profile \r'
        python3 general/performance_profile.py $RESFILES --outdir ${PLOTDIR}/performance_profile/ $PERFPROFARGS > /dev/null 2>&1; fi
        if test $TIME = "true"; then echo -ne '|░░░░░░░░            |  (40%)  Time Plot: Parsing      \r'
        for j in $OUTFILES; do python3 general/time.py $j --outdir $PKLDIR --save $TIMEARGS > /dev/null 2>&1; done; fi
        if test $TIME = "true"; then echo -ne '|░░░░░░░░░░░░        |  (60%)  Time Plot: Plotting     \r'
        python3 general/time.py $PKLDIR/check.*.pkl --outdir $PLOTDIR/timedist/ --compare $TIMEARGS > /dev/null 2>&1; fi
        if test $DETECTION = "true"; then echo -ne '|░░░░░░░░░░░░        |  (70%)  Detection Plot          \r'
        python3 detection/plotter_detection.py $OUTFILES --outdir $PLOTDIR/detection/ $DETECTIONARGS > /dev/null 2>&1; fi
        # The following scripts require GCG to be compiled with STATISTICS=true
        if test $LAST_STATISTICS = "true"; then
            if test $BOUNDS = "true"; then   echo -ne '|░░░░░░░░░░░░░░░░    |  (80%)  Bounds Plot         \r'
            #echo "calling python3 bounds/plotter_bounds.py $OUTFILES --outdir ${PLOTDIR}/bounds $BOUNDSARGS"
            python3 bounds/plotter_bounds.py $OUTFILES --outdir ${PLOTDIR}/bounds $BOUNDSARGS > /dev/null 2>&1; fi
        fi
        cd ../check
        echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done generating visualizations     \r'
        echo -ne '\n'
    fi
}

function printTeX_makefile(){
cat > comparisonreport.make << EndOfMessage

# latexmk automatically manages the .tex files
comparisonreport.pdf: ${REPORTFILE}
	@echo ------------
	@echo
	@echo Compiling tex code. This may take a while.
	@echo
	@echo ------------
	@latexmk -f -pdf -pdflatex="pdflatex -interaction=batchmode -shell-escape" -use-make ${REPORTFILE}
	@make -f comparisonreport.make clean

clean:
	@latexmk -c
	@rm -f report_*figure*.*
	@rm -f *.auxlock
	@rm -f *figure*.md5
	@rm -f *figure*.log
	@rm -f *figure*.dpth

cleanall:
	@latexmk -C
	@make -f comparisonreport.make clean
EndOfMessage
}

function printTeX_readme(){
cat > comparisonreport.README << EndOfMessage
README: How to create a PDF file from the .tex file(s) using the ${TSTNAME}.make file.
Note: The package pdflatex is required.

Use the command
    'make -f comparisonreport.make'
to compile.
Depending on the size of your problem that may take some time.
Please do not delete any new files that might be generated during the compile process.
All access files will be deleted automatically once the compilation is complete.

Clean options:
    'make -f comparisonreport.make clean' clears all present intermediate files (if any exist)
    'make -f comparisonreport.make cleanall' clears all generated files INCLUDING .pdf
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
% * Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       *
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
\usepackage{float}
\usepackage[pdfpagelabels]{hyperref}
\usepackage{fancybox}
\usepackage{longtable}
\usepackage{array,booktabs,color}
\setcounter{secnumdepth}{4}
\usepackage{graphicx}
\usepackage{grffile}

\begin{document}

\begin{titlepage}
  \begin{center}
  \thispagestyle{empty}
  {\Huge Comparison Report} \\\\ \today

\vspace{2cm}

\end{center}
\begin{tabular}{{lp{10cm}}}
    \multicolumn{2}{l}{\large\textbf{Test Run Characteristics Overview}} \\\\
EndOfMessage
}

function printTeX_TOC(){
cat >> ${REPORTFILE} << EndOfMessage
\thispagestyle{empty}
\tableofcontents
\newpage
EndOfMessage
}

function printTeX_frontpage_runs(){
  echo "\textbf{Number of Runs}:& ${NRUNS} \\\\" >> ${REPORTFILE}
  i=1
  for o in $OUTFILES; do
    o=$(basename $o)
    #outfiles[$i]=${o%.out}
    echo "\textbf{Run $i}:& ${o%.out} \\\\" | sed 's/\_/\\\_/g' >> ${REPORTFILE}
    ((i=i+1))
  done

cat >> ${REPORTFILE} << EndOfMessage
\end{tabular}
\end{titlepage}
\newpage
EndOfMessage
}

function printTeX_frontpage_table(){
# end table
echo "\end{tabular}" >> $REPORTFILE
# write table on front page
cat < ${PLOTDIR}/*.tex | head -n -1 | tail -n +5 | sed -e "s/{tabular\*}{\\\columnwidth}/{longtable}/g" -e "s/{tabular\*}/{longtable}/g" >> ${REPORTFILE}
# end front page
cat >> ${REPORTFILE} << EndOfMessage
\end{titlepage}
\newpage
EndOfMessage
# write TOC
printTeX_TOC;
}

function printTeX_secondpage_table(){
# write front page with runs
printTeX_frontpage_runs;
# write table on second page
cat < ${PLOTDIR}/*.tex | head -n -1 | tail -n +5 | sed -e "s/{tabular\*}{\\\columnwidth}/{longtable}/g" -e "s/{tabular\*}/{longtable}/g" >> ${REPORTFILE}
# end second page
echo "\\newpage" >> $REPORTFILE
# write TOC
printTeX_TOC;
}

# define descriptions for each visualization (deactivate automatic line breaks in your IDE for improved readibility)
declare -A desc_comp
desc_comp["performance_profile"]="The best performing run is a line with steepness \$m=0\$ and \$y\$-intercept \$f(0)=1\$. For all other runs, the integral over the shown \$x\$-interval is lower or equal to the one of the best performing run. The plot can be interpreted as follows: For each point \$(x,y)\$, \$x\$ represents the performance in multiples of the duration the best performing run needed and \$y\$ represents the likelihood for the setting to perform as fast as the best performing run times the multiples. If the underlying function of one setting has first-order stochastic dominance over another, it can be considered to be faster for this test set."
desc_comp["averagesolvetime"]="Average solving time."
desc_comp["failcomparison"]="Comparison of Fails."
desc_comp["runtimecomparison"]="Version to version runtime comparison."
desc_comp["runtimes"]="Runtimes per test run."
desc_comp["timecomparisonperstatus"]="Time comparison per status."
desc_comp["compare"]="Time Comparison Bar Chart."
desc_comp["comparepie"]="Time Comparison Pie Chart."
desc_comp["bounds_time"]="The top subplot shows the development of the primal and dual bounds in the RMP during the pricing in the root node as given by the table ``root bounds'' printed by GCG. Every change represents a pricing iteration and the resulting changes to the bounds. The bounds are complemented by a gap plot. The other two subplots illustrate the point in time in the pricing at which the columns that are finally in the basis are generated."
desc_comp["detection-times"]="This plot visualizes what fraction of instances against the time used for the whole detection process, including classification and score computation, on a logarithmic scale."
desc_comp["detection-decomps"]="The above plot shows the number of decompositions that were found for what fraction of instances in the given test set on a logarithmic scale. Only those decompositions are shown that have a score that is strictly greater than \$0\$."
desc_comp["detection-quality"]="This plot illustrates the detection goodness of the test set (according to the ``max white score''). If for all instances in the test set, the best (``whitest'') decomposition was completely white, this plot would show a line with \$y=1\$ for all \$x\$."
desc_comp["detection-quality_custom"]="Similarly to the previous plot, this one shows the detection goodness of the whole test set, but for a specific detector, the ``$CUSTOMDETECTOR'' detector. For this specific detector, GCG outputs the scores separately in the \texttt{detectionstatistics} test mode. With this plot, together with the previous one, one can compare the performance of the Set Partitioning Master detector with the overall performance."
desc_comp["detection-nBlocksOfBest"]="This plot shows how many blocks are used in the (according to the max white score) best decomposition, allowing to make statements about the sensibility of different numbers of blocks."
desc_comp["detection-classification_classes_custom"]="With this visualization, one can see how many classes a classifier determined for the variables. It can be chosen for which classifier to plot this, here ``$CUSTOMCLASSIFIER'' was chosen."

function printTeX_comparison_information(){
# start generating aggregated
# types and subtypes of visus:
# - performance profile (1 total)
# - general plots (5 total)
# - time distribution (2 total)
# - detection plots (6 total)
# - bounds plot (1 total)

declare -A titles
titles["performance_profile"]="Performance Profile"
titles["general"]="General Visualizations"
titles["timedist"]="Time Distribution Plot"
titles["detection"]="Detection Visualizations"
titles["bounds"]="Bounds Development (Root Node)"

visus=''
if [[ $PERFPROF = "true" ]]; then visus="$visus performance_profile"; fi
if [[ $GENERAL = "true" ]]; then visus="$visus general"; fi
if [[ $TIME = "true" ]]; then visus="$visus timedist"; fi
if [[ $DETECTION = "true" ]]; then visus="$visus detection"; fi
if [[ ($BOUNDS = "true") && ($LAST_STATISTICS = "true") ]]; then visus="$visus bounds"; fi

#echo "\newpage\section{Runtime Visualizations}" >> ${REPORTFILE}
for p in $visus
do
  echo "\section{${titles["${p}"]}}" >> ${REPORTFILE}
  #echo "\addcontentsline{toc}{section}{${titles["${p}"]}}" >> ${REPORTFILE}
  for v in $(ls plots/$p/*.pdf | sed "/bounds\.time/d") # subtype of visu
  do
    if [[ $p = "bounds" ]]; then echo "\subsection{Instance: $(echo $v | rev | cut -d "/" -f1 | rev | cut -d "." -f1 | sed "s/\_/\\\_/g")}" >> ${REPORTFILE}; fi
    echo "\vspace*{\fill}"    >> ${REPORTFILE}
    echo "\begin{figure}[H]"  >> ${REPORTFILE}
    echo "\makebox[\textwidth][c]{\includegraphics[width=1.2\textwidth]{$v}}" >> ${REPORTFILE}
    echo "  \caption{${desc_comp[$(echo $v | rev | cut -d "/" -f1 | cut -d "." -f2 | rev)]} \\\\ Visualization Path: \texttt{$(echo $v | sed 's/\_/\\\_/g')}}" >> ${REPORTFILE}
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
  printTeX_secondpage_table
  echo -ne '|░░░░░░              |  (30%)  Generating TeX Report Code... \r'
  if [ $NRUNS -gt 1 ]; then
    printTeX_comparison_information > /dev/null 2>&1
  fi
  printTeX_tail;
  echo -ne '|░░░░░░░░░░░░        |  (60%)  Compiling TeX Report Code...  \r'
  make -f comparisonreport.make > /dev/null 2>&1
  echo -ne '|░░░░░░░░░░░░░░░░░░░░|  (100%) Done! Opening Report...       \r'
  echo -ne '\n'
  xdg-open comparisonreport.pdf > /dev/null 2>&1
  if [ $? -eq 1 ]; then echo "Opening PDF.";
  else echo "Could not open PDF file."; fi
}

function main(){
  echoBold "Starting Comparison Report Generation...\n"
  extractScriptSettings;
  prepareVisualizationGeneration;
  echoBold "Generating visualizations."
  generateVisualizations;
  echoBold "Generating PDF report file.\n"
  generateReport;
}

main;
