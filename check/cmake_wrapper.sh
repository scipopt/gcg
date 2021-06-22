#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2021 Operations Research, RWTH Aachen University       *
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
# @author Erik Muehmer
#

if [[ -z $1 ]]
then
    exit 1;
fi

if [[ $1 = "cluster" ]]
then
    if [[ -z ${TSTNAME} ]]; then TSTNAME="short"; fi
    if [[ -z ${BINNAME} ]]; then BINNAME=${GCG_BINARY}; fi
    if [[ -z ${SETNAME} ]]; then SETNAME="default"; fi
    if [[ -z ${MSETNAME} ]]; then MSETNAME="default"; fi
    if [[ -z ${BINID} ]]; then BINID="1"; fi
    if [[ -z ${TIMELIMIT} ]]; then TIMELIMIT="3600"; fi
    if [[ -z ${NODELIMIT} ]]; then NODELIMIT="-1"; fi
    if [[ -z ${MEMLIMIT} ]]; then MEMLIMIT="16000"; fi
    if [[ -z ${THREADS} ]]; then THREADS="0"; fi
    if [[ -z ${FEASTOL} ]]; then FEASTOL="default"; fi
    if [[ -z ${LPS} ]]; then LPS="spx"; fi
    if [[ -z ${DISPFREQ} ]]; then DISPFREQ="100"; fi
    if [[ -z ${CONTINUE} ]]; then CONTINUE="false"; fi
    if [[ -z ${QUEUETYPE} ]]; then QUEUETYPE="condor"; fi
    if [[ -z ${QUEUE} ]]; then QUEUE="normal"; fi
    if [[ -z ${PPN} ]]; then PPN="1"; fi
    if [[ -z ${CLIENTTMPDIR} ]]; then CLIENTTMPDIR="/tmp/"; fi
    if [[ -z ${NOWAITCLUSTER} ]]; then NOWAITCLUSTER="1"; fi
    if [[ -z ${EXCLUSIVE} ]]; then EXCLUSIVE="notneeded"; fi
    if [[ -z ${PERMUTE} ]]; then PERMUTE="0"; fi
    if [[ -z ${MODE} ]]; then MODE="readdec"; fi
    if [[ -z ${STATISTICS} ]]; then STATISTICS="false"; fi
    if [[ -z ${PROJECT} ]]; then PROJECT="none"; fi
    if [[ -z ${SSH_USER} ]]; then SSH_USER=$USER; fi

    # print parameters
    echo "Parameters:"
    echo "TSTNAME       = $TSTNAME"
    echo "BINNAME       = $BINNAME"
    echo "SETNAME       = $SETNAME"
    echo "MSETNAME      = $MSETNAME"
    echo "BINID         = $BINID"
    echo "TIMELIMIT     = $TIMELIMIT"
    echo "NODELIMIT     = $NODELIMIT"
    echo "MEMLIMIT      = $MEMLIMIT"
    echo "THREADS       = $THREADS"
    echo "FEASTOL       = $FEASTOL"
    echo "LPS           = $LPS"
    echo "DISPFREQ      = $DISPFREQ"
    echo "CONTINUE      = $CONTINUE"
    echo "QUEUETYPE     = $QUEUETYPE"
    echo "QUEUE         = $QUEUE"
    echo "PPN           = $PPN"
    echo "CLIENTTMPDIR  = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER = $NOWAITCLUSTER"
    echo "EXCLUSIVE     = $EXCLUSIVE"
    echo "PERMUTE       = $PERMUTE"
    echo "MODE          = $MODE"
    echo "STATISTICS    = $STATISTICS"
    echo "PROJECT       = $PROJECT"
    echo "SSH_USER      = $SSH_USER"

    if [[ $QUEUETYPE = "condor" ]]
    then
        echo "Connecting to clustor.or.rwth-aachen.de to submit jobs ..."
        ssh -t "$SSH_USER"@clustor.or.rwth-aachen.de "cd $PWD ; ./check_cluster.sh \"${TSTNAME}\" \"${BINNAME}\" \"${SETNAME}\" \"${MSETNAME}\" \"${BINID}\" \"${TIMELIMIT}\" \"${NODELIMIT}\" \"${MEMLIMIT}\" \"${THREADS}\" \"${FEASTOL}\" \"${LPS}\" \"${DISPFREQ}\" \"${CONTINUE}\" \"${QUEUETYPE}\" \"${QUEUE}\" \"${PPN}\" \"${CLIENTTMPDIR}\" \"${NOWAITCLUSTER}\" \"${EXCLUSIVE}\" \"${PERMUTE}\" \"${MODE}\" \"${STATISTICS}\" \"${PROJECT}\""
     else
        . ./check_cluster.sh "${TSTNAME}" "${BINNAME}" "${SETNAME}" "${MSETNAME}" "${BINID}" "${TIMELIMIT}" "${NODELIMIT}" "${MEMLIMIT}" "${THREADS}" "${FEASTOL}" "${LPS}" "${DISPFREQ}" "${CONTINUE}" "${QUEUETYPE}" "${QUEUE}" "${PPN}" "${CLIENTTMPDIR}" "${NOWAITCLUSTER}" "${EXCLUSIVE}" "${PERMUTE}" "${MODE}" "${STATISTICS}" "${PROJECT}"
     fi
fi
