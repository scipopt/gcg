#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
# @author Erik Muehmer
#

if [[ -z $1 ]]
then
    exit 1;
fi

if [[ $1 = "visu" ]]
then
    if [[ -z ${VISUSETTINGS} ]]; then VISUSETTINGS="none"; fi
    if [[ -z ${BINID} ]]; then BINID=unknown; fi
    if [[ -z ${VERSION} ]]; then VERSION=unkown; fi
    if [[ -z ${MODE} ]]; then MODE=unkown; fi
    if [[ -z ${LPS} ]]; then LPS=unkown; fi
    if [[ -z ${THREADS} ]]; then THREADS=unkown; fi
    if [[ -z ${FEASTOL} ]]; then FEASTOL=unkown; fi
    if [[ -z ${LAST_STATISTICS} ]]; then LAST_STATISTICS="true"; fi

    ./writeTestsetReport.sh $VISUSETTINGS $BINID $VERSION $MODE $LPS $THREADS $FEASTOL $LAST_STATISTICS
fi

if [[ $1 = "cluster" ]]
then
    if [[ -z ${TEST} ]]; then TEST="short"; fi
    if [[ -z ${BINNAME} ]]; then BINNAME=${GCG_BINARY}; fi
    if [[ -z ${SETTINGS} ]]; then SETTINGS="default"; fi
    if [[ -z ${MASTERSETTINGS} ]]; then MASTERSETTINGS="default"; fi
    if [[ -z ${BINID} ]]; then BINID="1"; fi
    if [[ -z ${TIME} ]]; then TIME="3600"; fi
    if [[ -z ${NODES} ]]; then NODES="-1"; fi
    if [[ -z ${MEM} ]]; then MEM="16000"; fi
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
    echo "TEST              = $TEST"
    echo "BINNAME           = $BINNAME"
    echo "SETTINGS          = $SETTINGS"
    echo "MASTERSETTINGS    = $MASTERSETTINGS"
    echo "BINID             = $BINID"
    echo "TIME              = $TIME"
    echo "NODES             = $NODES"
    echo "MEM               = $MEM"
    echo "THREADS           = $THREADS"
    echo "FEASTOL           = $FEASTOL"
    echo "LPS               = $LPS"
    echo "DISPFREQ          = $DISPFREQ"
    echo "CONTINUE          = $CONTINUE"
    echo "QUEUETYPE         = $QUEUETYPE"
    echo "QUEUE             = $QUEUE"
    echo "PPN               = $PPN"
    echo "CLIENTTMPDIR      = $CLIENTTMPDIR"
    echo "NOWAITCLUSTER     = $NOWAITCLUSTER"
    echo "EXCLUSIVE         = $EXCLUSIVE"
    echo "PERMUTE           = $PERMUTE"
    echo "MODE              = $MODE"
    echo "STATISTICS        = $STATISTICS"
    echo "PROJECT           = $PROJECT"
    echo "SSH_USER          = $SSH_USER"

    if [[ $QUEUETYPE = "condor" ]]
    then
        echo "Connecting to clustor.or.rwth-aachen.de to submit jobs ..."
        ssh -t "$SSH_USER"@clustor.or.rwth-aachen.de "cd $PWD ; ./check_cluster.sh \"${TEST}\" \"${BINNAME}\" \"${SETTINGS}\" \"${MASTERSETTINGS}\" \"${BINID}\" \"${TIME}\" \"${NODES}\" \"${MEM}\" \"${THREADS}\" \"${FEASTOL}\" \"${LPS}\" \"${DISPFREQ}\" \"${CONTINUE}\" \"${QUEUETYPE}\" \"${QUEUE}\" \"${PPN}\" \"${CLIENTTMPDIR}\" \"${NOWAITCLUSTER}\" \"${EXCLUSIVE}\" \"${PERMUTE}\" \"${MODE}\" \"${STATISTICS}\" \"${PROJECT}\""
     else
        . ./check_cluster.sh "${TEST}" "${BINNAME}" "${SETTINGS}" "${MASTERSETTINGS}" "${BINID}" "${TIME}" "${NODES}" "${MEM}" "${THREADS}" "${FEASTOL}" "${LPS}" "${DISPFREQ}" "${CONTINUE}" "${QUEUETYPE}" "${QUEUE}" "${PPN}" "${CLIENTTMPDIR}" "${NOWAITCLUSTER}" "${EXCLUSIVE}" "${PERMUTE}" "${MODE}" "${STATISTICS}" "${PROJECT}"
     fi
fi

if [[ $1 = "test" ]]
then
    if [[ -z ${TEST} ]]; then TEST="short"; fi
    if [[ -z ${BINNAME} ]]; then BINNAME=${GCG_BINARY}; fi
    if [[ -z ${SETTINGS} ]]; then SETTINGS="default"; fi
    if [[ -z ${MASTERSETTINGS} ]]; then MASTERSETTINGS="default"; fi
    if [[ -z ${BINID} ]]; then BINID="1"; fi
    if [[ -z ${TIME} ]]; then TIME="3600"; fi
    if [[ -z ${NODES} ]]; then NODES="-1"; fi
    if [[ -z ${MEM} ]]; then MEM="16000"; fi
    if [[ -z ${THREADS} ]]; then THREADS="0"; fi
    if [[ -z ${FEASTOL} ]]; then FEASTOL="default"; fi
    if [[ -z ${DISPFREQ} ]]; then DISPFREQ="100"; fi
    if [[ -z ${CONTINUE} ]]; then CONTINUE="false"; fi
    if [[ -z ${LOCK} ]]; then LOCK="false"; fi
    if [[ -z ${VERSION} ]]; then VERSION="unknown"; fi
    if [[ -z ${LPS} ]]; then LPS="spx"; fi
    if [[ -z ${VALGRIND} ]]; then VALGRIND="false"; fi
    if [[ -z ${MODE} ]]; then MODE="readdec"; fi
    if [[ -z ${SETCUTOFF} ]]; then SETCUTOFF="false"; fi
    if [[ -z ${STATISTICS} ]]; then STATISTICS="false"; fi
    if [[ -z ${SHARED} ]]; then SHARED="false"; fi
    if [[ -z ${VISU} ]]; then VISU="false"; fi
    if [[ -z ${LAST_STATISTICS} ]]; then LAST_STATISTICS="false"; fi
    if [[ -z ${VISUSETTINGS} ]]; then VISUSETTINGS="none"; fi
    if [[ -z ${DETECTIONSTATISTICS} ]]; then DETECTIONSTATISTICS="false"; fi

    # print parameters
    echo "Parameters:"
    echo "TEST              = $TEST"
    echo "BINNAME           = $BINNAME"
    echo "SETTINGS          = $SETTINGS"
    echo "MASTERSETTINGS    = $MASTERSETTINGS"
    echo "BINID             = $BINID"
    echo "TIME              = $TIME"
    echo "NODES             = $NODES"
    echo "MEM               = $MEM"
    echo "THREADS           = $THREADS"
    echo "FEASTOL           = $FEASTOL"
    echo "DISPFREQ          = $DISPFREQ"
    echo "CONTINUE          = $CONTINUE"
    echo "LOCK              = $LOCK"
    echo "VERSION           = $VERSION"
    echo "LPS               = $LPS"
    echo "VALGRIND          = $VALGRIND"
    echo "MODE              = $MODE"
    echo "SETCUTOFF         = $SETCUTOFF"
    echo "STATISTICS        = $STATISTICS"
    echo "SHARED            = $SHARED"
    echo "VISU              = $VISU"
    echo "LAST_STATISTICS   = $LAST_STATISTICS"
    echo "VISUSETTINGS      = $VISUSETTINGS"
    echo "DETECTIONSTATISTICS = $DETECTIONSTATISTICS"

    . ./check.sh "${TEST}" "${BINNAME}" "${SETTINGS}" "${MASTERSETTINGS}" "${BINID}" "${TIME}" "${NODES}" "${MEM}" "${THREADS}" "${FEASTOL}" "${DISPFREQ}" "${CONTINUE}" "${LOCK}" "${VERSION}" "${LPS}" "${VALGRIND}" "${MODE}" "${SETCUTOFF}" "${STATISTICS}" "${SHARED}" "${VISU}" "${LAST_STATISTICS}" "${VISUSETTINGS}" "${DETECTIONSTATISTICS}"
fi
