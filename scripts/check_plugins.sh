#!/usr/bin/env bash
#
# Check which SCIP plugins are included in GCG
# and print those which are not yet included.
#
# Invocation: scripts/check_plugins.sh INCLUDEFILE
#
# where INCLUDEFILE is a .c source code file
# defining which plugins are included in GCG
# (e.g. gcgplugins.c or masterplugins.c)
#

INCLUDEFILE=$1

SCIPINCLUDE=lib/scip-git/src/scip/scipdefplugins.c

# List of plugins which should be ignored
EXCLUDE="
SCIPincludeConshdlrNonlinear
SCIPincludeConshdlrQuadratic
SCIPincludeDispDefault
SCIPincludeHeurNlpdiving
SCIPincludeHeurTrivial
SCIPincludeHeurTrysol
SCIPincludeHeurUndercover
SCIPincludeNodeselBfs
SCIPincludeNodeselDfs
SCIPincludeNodeselEstimate
SCIPincludeNodeselHybridestim
SCIPincludeNodeselRestartdfs
SCIPincludePresolComponents
SCIPincludePresolConvertinttobin
SCIPincludePropObbt
"


PLUGINS=`grep -e '^ *SCIP_CALL( *SCIPinclude[A-Za-z0-9]*(scip) *)' $SCIPINCLUDE |
    sed 's/.*\(SCIPinclude[A-Za-z0-9]*\).*/\1/'`

for PLUGIN in $PLUGINS
do
    if [ `grep -c $PLUGIN $INCLUDEFILE` -eq 0 ] &&
        [ `echo $EXCLUDE | grep -c $PLUGIN` -eq 0 ]
    then
        echo $PLUGIN
    fi
done
