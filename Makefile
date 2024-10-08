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

#@file    Makefile
#@brief   Makefile for generic column generation code using SCIP as a callable library
#@author  Gerald Gamrath
#@author  Martin Bergner
#@author  Christian Puchert
#@author  Matthias Walter

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------
VERSION         :=	3.7.0
GCGGITHASH	=
SCIPDIR         =   $(CURDIR)/lib/scip

#-----------------------------------------------------------------------------
# necessary information
#-----------------------------------------------------------------------------
LIBDIR          =	lib
DIRECTORIES     =	$(LIBDIR) $(LIBDIR)/shared $(LIBDIR)/include $(LIBDIR)/static $(LIBOBJDIR) $(LIBOBJSUBDIRS)
SOFTLINKS	=
MAKESOFTLINKS	=	true

SHELL		= 	bash
READ		=	read -e
LN_s		= 	ln -s
GCGDIR		=	$(realpath .)
TIME		=	3600
DIP			=	dip
MASTERSETTINGS	=	default
VISUSETTINGS 	= 	none
DATADIR 		= 	none
OBJDIR		= obj

VALGRIND	=	false
MODE		=	readdec
VISU 		= 	false
DETECTIONSTATISTICS = 	false
STATISTICS  =  	false
PROJECT		=	none
GTEST		=	false
PARASCIP	= 	true
SYM      	=   snauty
CLIQUER     =   false
HMETIS      =   false
OPENMP      =   false
GSL         =   false
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings
LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(BLISS).$(CLIQUER)

# overriding SCIP PARASCIP setting if compiled with OPENMP
ifeq ($(OPENMP),true)
override PARASCIP=true
MAKEOVERRIDES += PARASCIP=true
endif

# overriding SCIP LPS setting if compiled with CPLEXSOLVER
ifeq ($(CPLEXSOLVER),true)
override LPS =  cpx
MAKEOVERRIDES += LPS=cpx
endif

# overriding sym settings
ifneq (,$(filter $(SYM),bliss sbliss))
override BLISS=true
override NAUTY=false
MAKEOVERRIDES += BLISS=true
MAKEOVERRIDES += NAUTY=false
endif

ifneq (,$(filter $(SYM),nauty snauty))
override BLISS=false
override NAUTY=true
MAKEOVERRIDES += BLISS=false
MAKEOVERRIDES += NAUTY=true
endif

ifeq ($(SYM),none)
override BLISS=false
override NAUTY=false
MAKEOVERRIDES += BLISS=false
MAKEOVERRIDES += NAUTY=false
endif

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP (need to do this twice, once to
# find the correct binary, then, after getting the correct flags from the
# binary (which is necessary since the ZIMPL flags differ from the default
# if compiled with the SCIP Optsuite instead of SCIP), we need to set the
# compile flags, e.g., for the ZIMPL library, which is again done in make.project.
# If the include fails (silently), then the code below the inclusion is not
# supposed to be run.
#-----------------------------------------------------------------------------
-include $(SCIPDIR)/make/make.project
ifneq ($(wildcard $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION)), )
SCIPVERSION       :=$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override ARCH     :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ARCH=\([^@]*\).*/\1/')
override EXPRINT  :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* EXPRINT=\([^@]*\).*/\1/')
override GAMS     :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GAMS=\([^@]*\).*/\1/')
override GMP      :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM      :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT    :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPSCHECK :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT   :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override PARASCIP :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* PARASCIP=\([^@]*\).*/\1/')
override READLINE :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override ZIMPL    :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB     :=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
include $(SCIPDIR)/make/make.project
endif

#-----------------------------------------------------------------------------
# BLISS
#-----------------------------------------------------------------------------

ifeq ($(BLISS),true)
FLAGS		+=	-DWITH_BLISS
LDFLAGS		+= 	-lbliss
ifeq ($(COMP),gnu)
FLAGS		+=	-isystem$(LIBDIR)/include
else
FLAGS		+=	-I$(LIBDIR)/include
endif
SOFTLINKS	+=	$(LIBDIR)/include/bliss
SOFTLINKS	+=	$(LIBDIR)/static/libbliss.$(STATICLIBEXT)
LINKMSG		+=	"bliss graph isomorphism framework (disable by compiling with \"make BLISS=false\"):\n"
LINKMSG		+=	" -> bliss is the path to the bliss include files, e.g., \"bliss-0.72\"\n"
LINKMSG		+=	" -> \"libbliss.$(STATICLIBEXT)\" is the path to the bliss library, e.g., \"bliss-0.72/libbliss.$(STATICLIBEXT)\"\n"
endif

#-----------------------------------------------------------------------------
# Cliquer
#-----------------------------------------------------------------------------

ifeq ($(CLIQUER),true)
FLAGS		+=	-DWITH_CLIQUER
LDFLAGS		+= 	-lcliquer
ifeq ($(COMP),gnu)
FLAGS		+=	-isystem$(LIBDIR)/include
else
FLAGS		+=	-I$(LIBDIR)/include
endif
SOFTLINKS	+=	$(LIBDIR)/include/cliquer
SOFTLINKS	+=	$(LIBDIR)/static/libcliquer.$(STATICLIBEXT)
LINKMSG		+=	"cliquer library (disable by compiling with \"make CLIQUER=false\"):\n"
LINKMSG		+=	" -> cliquer is the path to the cliquer include files, e.g., \"cliquer-1.21\"\n"
LINKMSG		+=	" -> \"libcliquer.$(STATICLIBEXT)\" is the path to the cliquer library, e.g., \"cliquerinc/libcliquer.$(STATICLIBEXT)\"\n"
endif

#-----------------------------------------------------------------------------
# hmetis
#-----------------------------------------------------------------------------

ifeq ($(HMETIS),true)
FLAGS		+=	-DWITH_HMETIS
endif

#-----------------------------------------------------------------------------
# GSL
#-----------------------------------------------------------------------------
ifeq ($(GSL),true)
LDFLAGS                +=      -lgsl -lgslcblas -lm
FLAGS		           +=	   -DWITH_GSL
endif

#-----------------------------------------------------------------------------
# CPLEX pricing solver
#-----------------------------------------------------------------------------

ifeq ($(CPLEXSOLVER),true)
FLAGS		+=	-DWITH_CPLEXSOLVER -I$(SCIPDIR)/lib/include/cpxinc
endif

#-----------------------------------------------------------------------------
# GCG statistics
#-----------------------------------------------------------------------------

ifeq ($(STATISTICS),true)
FLAGS		+=	-DSCIP_STATISTIC
endif

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gcg

LIBOBJ = \
			gcg/branch_empty.o \
			gcg/branch_generic.o \
			gcg/benders_gcg.o \
			gcg/branch_orig.o \
			gcg/branch_relpsprob.o \
			gcg/branch_ryanfoster.o \
			gcg/branch_bpstrong.o \
			gcg/class_conspartition.o \
			gcg/class_indexpartition.o \
			gcg/class_pricingcontroller.o \
			gcg/class_pricingtype.o \
			gcg/class_partialdecomp.o \
			gcg/class_detprobdata.o \
			gcg/class_stabilization.o \
			gcg/class_varpartition.o \
			gcg/clscons_miplibconstypes.o \
			gcg/clscons_nnonzeros.o \
			gcg/clscons_consnamenonumbers.o \
			gcg/clscons_consnamelevenshtein.o \
			gcg/clscons_scipconstypes.o \
			gcg/clscons_gamsdomain.o \
			gcg/clscons_gamssymbol.o \
			gcg/clscons.o \
			gcg/clsvar_gamsdomain.o \
			gcg/clsvar_gamssymbol.o \
			gcg/clsvar_objvalues.o \
			gcg/clsvar_scipvartypes.o \
			gcg/clsvar_objvaluesigns.o \
			gcg/clsvar.o \
			gcg/colpool.o \
			gcg/cons_decomp.o \
			gcg/cons_integralorig.o \
			gcg/cons_masterbranch.o \
			gcg/cons_origbranch.o \
			gcg/dec_compgreedily.o \
			gcg/dec_connected_noNewLinkingVars.o \
			gcg/dec_connectedbase.o \
			gcg/dec_consclass.o \
			gcg/dec_constype.o \
			gcg/dec_dbscan.o \
			gcg/dec_densemasterconss.o \
			gcg/dec_generalmastersetcover.o \
			gcg/dec_generalmastersetpack.o \
			gcg/dec_generalmastersetpart.o \
			gcg/dec_hcgpartition.o \
			gcg/dec_hrcgpartition.o \
			gcg/dec_hrgpartition.o \
			gcg/dec_mastersetcover.o \
			gcg/dec_mastersetpack.o \
			gcg/dec_mastersetpart.o \
			gcg/dec_mst.o \
			gcg/dec_neighborhoodmaster.o \
			gcg/dec_postprocess.o \
			gcg/dec_staircase_lsp.o \
			gcg/dec_stairheur.o \
			gcg/dec_varclass.o \
			gcg/decomp.o \
			gcg/dialog_explore.o \
			gcg/dialog_gcg.o \
			gcg/dialog_graph.o \
			gcg/dialog_master.o \
			gcg/disp_gcg.o \
			gcg/disp_master.o \
			gcg/event_bestsol.o \
			gcg/event_display.o \
			gcg/event_mastersol.o \
			gcg/event_relaxsol.o \
			gcg/event_solvingstats.o \
			gcg/gcgcol.o \
			gcg/gcg_general.o \
			gcg/gcggithash.o \
			gcg/gcgheur.o \
			gcg/gcgsepa.o \
			gcg/gcgplugins.o \
			gcg/gcgpqueue.o \
			gcg/gcgsort.o \
			gcg/gcgvar.o \
			graph/graph_gcg.o \
			graph/graph_tclique.o \
			graph/inst.o \
			graph/weights.o \
			gcg/heur_gcgcoefdiving.o \
			gcg/heur_gcgdins.o \
			gcg/heur_gcgfeaspump.o \
			gcg/heur_gcgfracdiving.o \
			gcg/heur_gcgguideddiving.o \
			gcg/heur_gcglinesdiving.o \
			gcg/heur_gcgpscostdiving.o \
			gcg/heur_gcgrens.o \
			gcg/heur_gcgrins.o \
			gcg/heur_gcgrounding.o \
			gcg/heur_gcgshifting.o \
			gcg/heur_gcgsimplerounding.o \
			gcg/heur_gcgveclendiving.o \
			gcg/heur_gcgzirounding.o \
			gcg/heur_greedycolsel.o \
			gcg/heur_mastercoefdiving.o \
			gcg/heur_masterdiving.o \
			gcg/heur_masterfracdiving.o \
			gcg/heur_masterlinesdiving.o \
			gcg/heur_mastervecldiving.o \
			gcg/heur_origdiving.o \
			gcg/heur_relaxcolsel.o \
			gcg/heur_restmaster.o \
			gcg/heur_setcover.o \
			gcg/heur_xpcrossover.o \
			gcg/heur_xprins.o \
			gcg/masterplugins.o \
			gcg/bendersplugins.o \
			gcg/misc.o \
			gcg/miscvisualization.o \
			gcg/nodesel_master.o \
			gcg/objdialog.o \
			gcg/params_visu.o \
			gcg/presol_roundbound.o \
			gcg/pricer_gcg.o \
			gcg/pricestore_gcg.o \
			gcg/pricingjob.o \
			gcg/pricingprob.o \
			gcg/reader_blk.o \
			gcg/reader_cls.o \
			gcg/reader_dec.o \
			gcg/reader_gp.o \
			gcg/reader_ref.o \
			gcg/reader_tex.o \
			gcg/relax_gcg.o \
			gcg/scip_misc.o \
			gcg/score_bender.o \
			gcg/score_border.o \
			gcg/score_classic.o \
			gcg/score_fawh.o \
			gcg/score_forswh.o \
			gcg/score_maxwhite.o \
			gcg/score_spfawh.o \
			gcg/score_spfwh.o \
			gcg/score_strong.o \
			gcg/score.o \
			gcg/sepa_basis.o \
			gcg/sepa_master.o \
			gcg/solver.o \
			gcg/solver_knapsack.o \
			gcg/solver_mip.o \
			gcg/stat.o \

ifneq ($(SYM),none)
LIBOBJ		+=	symmetry/automorphism.o \
			gcg/dec_isomorph.o
endif

ifeq ($(BLISS),true)
LIBOBJ		+=	symmetry/automorphism_bliss.o
endif

ifeq ($(NAUTY),true)
LIBOBJ		+=	symmetry/automorphism_nauty.o \
			nauty/nauty.o \
			nauty/nautil.o \
			nauty/nausparse.o \
			nauty/schreier.o \
			nauty/naurng.o
FLAGS		+=	-DWITH_NAUTY -I$(SCIPDIR)/src/nauty
endif

ifeq ($(CLIQUER),true)
LIBOBJ		+=	gcg/solver_cliquer.o
endif

ifeq ($(CPLEXSOLVER),true)
LIBOBJ		+=	gcg/solver_cplex.o
endif

MAINOBJ		=	main.o

MAINSRC		=	$(filter $(wildcard $(SRCDIR)/*.c),$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*.cpp),$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp)))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME)-$(VERSION).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

# GCG Library
LIBOBJDIR	=	$(OBJDIR)/lib
OBJSUBDIRS	= 	gcg graph symmetry
ifeq ($(NAUTY),true)
OBJSUBDIRS	+=	nauty
endif
LIBOBJSUBDIRS   =       $(addprefix $(LIBOBJDIR)/,$(OBJSUBDIRS))

GCGLIBSHORTNAME =	gcg
GCGLIBNAME	=	$(GCGLIBSHORTNAME)-$(VERSION)

GCGLIBOBJ	=	${LIBOBJ}
GCGLIB		=	$(GCGLIBNAME).$(BASE)
ifeq ($(SHARED),true)
GCGLIBFILE	=	$(LIBDIR)/shared/lib$(GCGLIB).$(LIBEXT)
GCGLIBLINK	=	$(LIBDIR)/shared/lib$(GCGLIBSHORTNAME).$(BASE).$(LIBEXT)
GCGLIBSHORTLINK = 	$(LIBDIR)/shared/lib$(GCGLIBSHORTNAME).$(LIBEXT)
LDFLAGS		+=	$(LINKCXX_L)$(LIBDIR)/shared/
else
GCGLIBFILE	=	$(LIBDIR)/static/lib$(GCGLIB).$(LIBEXT)
GCGLIBLINK	=	$(LIBDIR)/static/lib$(GCGLIBSHORTNAME).$(BASE).$(LIBEXT)
GCGLIBSHORTLINK = 	$(LIBDIR)/static/lib$(GCGLIBSHORTNAME).$(LIBEXT)
LDFLAGS		+=	$(LINKCXX_L)$(LIBDIR)/static/
endif
GCGLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(GCGLIBOBJ))
GCGLIBSRC	=	$(filter $(wildcard $(SRCDIR)/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
GCGLIBSRC	+=	$(filter $(wildcard $(SRCDIR)/*/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
ifeq ($(NAUTY),true)
GCGLIBSRC	+=	$(filter $(wildcard $(SCIPDIR)/src/nauty/*.c),$(GCGLIBOBJ:.o=.c))
endif
GCGLIBDEP	=	$(SRCDIR)/depend.gcglib.$(OPT)

ALLSRC		=	$(MAINSRC) $(GCGLIBSRC)
SPLINT		=       splint
#SPLINTFLAGS	=	-UNDEBUG -UWITH_READLINE -UROUNDING_FE -UWITH_GMP -UWITH_ZLIB -preproc -formatcode +skip-sys-headers -weak +relaxtypes
SPLINTFLAGS	=	-UNDEBUG -UWITH_READLINE -UROUNDING_FE -UWITH_GMP -UWITH_ZLIB -which-lib -warn-posix-headers +skip-sys-headers -preproc -formatcode -weak \
			-redef +export-header +export-local +decl-undef +relaxtypes

GCGGITHASHFILE	= 	$(SRCDIR)/githash.c

#-----------------------------------------------------------------------------
# Flags
#-----------------------------------------------------------------------------

ifeq ($(OPENMP),true)
CFLAGS		+=	-fopenmp
LDFLAGS		+=	-fopenmp
CXXFLAGS	+=	-fopenmp
endif

ifeq ($(COMP),gnu)
CXXFLAGS	+=	-Wno-variadic-macros
endif

# WORKAROUND for missing DCXXFLAGS (C++ flags for dependency calls):
DCXXFLAGS=$(CXXFLAGS)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.PHONY: all
all:       	$(SCIPDIR)
		@-$(MAKE) libs
		@-$(MAKE) mainfiles

.PHONY: mainfiles
mainfiles:
		@-$(MAKE) $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

$(SCIPDIR)/make/make.project: |$(SCIPDIR)

-include make/local/make.targets

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK) ${GCGLIBFILE} ${GCGLIB} $(GCGLIBLINK) ${GCGLIBSHORTLINK} ${TESTSHORTLINK} ${GCGLIBOBJFILES} $(TESTOBJFILES) ${TESTFILE} ${TESTMAIN}
endif

.PHONY: libs
libs:		$(LINKSMARKERFILE) makegcglibfile $(GCGLIBLINK) $(GCGLIBSHORTLINK)

.PHONY: lint
lint:		$(ALLSRC)
		-rm -f lint.out
ifeq ($(FILES),)
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'
else
		$(SHELL) -ec  'for i in $(FILES); \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'
endif

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) $^ all

.PHONY: scip_clean
scip_clean:
		@$(MAKE) -C $(SCIPDIR) $^ clean

.PHONY: splint
splint:		$(ALLSRC)
		-rm -f splint.out
ifeq ($(FILES),)
		$(SHELL) -c '$(SPLINT) -I$(SRCDIR) -I/usr/include/linux $(FLAGS) $(SPLINTFLAGS)  $(filter %.c %.h,$^) &>> splint.out;'
else
		$(SHELL) -c '$(SPLINT) -I$(SRCDIR) -I/usr/include/linux $(FLAGS) $(SPLINTFLAGS) $(FILES %.c %.h,$^) &>> splint.out;'
endif

.PHONY: doc
doc:
		@cd doc; export BINDIR=$(PWD)/$(BINDIR); $(SHELL) builddoc.sh --mathjax;

.PHONY: $(MAINSHORTLINK)
$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

# include target to detect the current git hash
-include make/local/make.detectgithash

# this empty target is needed for the SCIP release versions
githash::   # do not remove the double-colon

${GCGGITHASHFILE}: githash

.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(MASTERSETTINGS) $(notdir $(BINDIR)/$(GCGLIBNAME).$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(MODE) $(SETCUTOFF)\
		$(STATISTICS) $(SHARED) $(VISU) $(LAST_STATISTICS) $(VISUSETTINGS) $(DETECTIONSTATISTICS);

.PHONY: visu
visu:
		cd check; \
		$(SHELL) ./writeTestsetReport.sh $(VISUSETTINGS) $(BINID) $(VERSION) $(MODE) $(LPS) $(THREADS) $(FEASTOL) $(LAST_STATISTICS) $(DATADIR);

.PHONY: eval
eval:
		cd check; \
		$(SHELL) ./eval.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(BINDIR)/$(GCGLIBNAME).$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND);


.PHONY: cleanbin
cleanbin:       $(BINDIR)
		@echo "-> remove binary $(MAINFILE)"
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: cleanlib
cleanlib:       $(LIBDIR)
		@echo "-> remove library $(GCGLIBFILE)"
		@-rm -f $(GCGLIBFILE) $(GCGLIBLINK) $(GCGLIBSHORTLINK)

.PHONY: cleantest
cleantest:
ifneq ($(OBJDIR),)
		@$(SHELL) -ec 'if test -d $(OBJDIR)/tests/; \
			then \
				echo "-> remove $(OBJDIR)/tests/"; \
				rm -f -f $(OBJDIR)/tests/*.o ; \
				cd $(OBJDIR) && rmdir tests ; \
			fi'
endif

.PHONY: clean
clean:          cleantest cleanlib cleanbin  $(DIRECTORIES)
		@-rm -f $(LASTSETTINGS)
ifneq ($(LIBOBJDIR),)
		@-(rm -f $(LIBOBJDIR)/*.o)
		@-(cd $(LIBOBJDIR) && rm -f */*.o && rmdir $(OBJSUBDIRS));
		@-rmdir $(LIBOBJDIR);
endif
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o);
		@-(rmdir $(OBJDIR));
endif


.PHONY: tags
tags:
		cd src/; rm -f TAGS; etags *.cpp *.c *.h ../$(SCIPDIR)/src/scip/*.c ../$(SCIPDIR)/src/scip/*.h ../$(SCIPDIR)/src/objscip/*.cpp ../$(SCIPDIR)/src/objscip/*.h;

.PHONY: depend
depend:		gcglibdepend testdepend | $(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(subst isystem,I,$(FLAGS)) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'
-include	$(MAINDEP)

.PHONY: gcglibdepend
gcglibdepend:
		$(SHELL) -ec '$(DCXX) $(DCXXFLAGS) $(subst isystem,I,$(FLAGS)) $(DFLAGS) $(GCGLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(GCGLIBDEP)'
-include	$(GCGLIBDEP)

.PHONY: testdepend
testdepend:: # do not remove double colon

tests:: #do not remove double colon

$(GCGLIBLINK):	$(GCGLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(GCGLIBFILE)) $(notdir $@)


$(MAINFILE):	$(SCIPLIBBASEFILE) $(LPILIBFILE) $(TPILIBFILE) $(MAINOBJFILES) $(GCGLIBFILE) | $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_l)$(GCGLIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_L)$(SCIPDIR)/lib/static $(LINKCXX_l)$(SCIPLIBBASE)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(TPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.c | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@

ifeq ($(NAUTY),true)
$(LIBOBJDIR)/nauty/%.o:	$(SCIPDIR)/src/nauty/%.c | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@
endif

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(OBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(OBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

.PHONY: makegcglibfile
makegcglibfile:  touchexternal $(GCGLIBFILE)

$(GCGLIBFILE):	$(LIBDIR) $(LIBDIR)/static $(LIBDIR)/shared $(GCGLIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(GCGLIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

.PHONY: $(GCGLIBSHORTLINK)
$(GCGLIBSHORTLINK):	$(GCGLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(GCGLIBFILE)) $(notdir $@)


$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		$(SOFTLINKS) | $(LIBDIR)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

$(DIRECTORIES):
		@echo
		@echo "- creating directory \"$@\""
		@-mkdir -p $@
-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal: | $(LIBOBJDIR)
ifneq ($(LAST_CPLEXSOLVER),$(CPLEXSOLVER))
		@-touch $(SRCDIR)/solver_cplex.c
endif
ifneq ($(LAST_STATISTICS),$(STATISTICS))
		@$(SHELL) -ec 'for file in $$(grep "SCIP_STATISTIC" ${SRCDIR} -rl | sort -u); do touch $$file; done'
endif

ifneq ($(LAST_BLISS),$(BLISS))
		@-touch $(SRCDIR)/dec_isomorph.cpp
		@-touch $(SRCDIR)/relax_gcg.c
		@-touch $(SRCDIR)/gcgplugins.c
endif
ifneq ($(LAST_CLIQUER),$(CLIQUER))
		@-touch $(SRCDIR)/solver_cliquer.c
		@-touch $(SRCDIR)/masterplugins.c
endif
ifneq ($(USRFLAGS),$(LAST_USRFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USROFLAGS),$(LAST_USROFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCFLAGS),$(LAST_USRCFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRLDFLAGS),$(LAST_USRLDFLAGS))
		@-touch -c $(GCGLIBOBJFILES) $(MAINOBJFILES)
endif
ifneq ($(USRARFLAGS),$(LAST_USRARFLAGS))
		@-touch -c $(GCGLIBOBJFILES) $(MAINOBJFILES)
endif
ifneq ($(OPENMP),$(LAST_OPENMP))
		@-touch $(ALLSRC)
endif
ifneq ($(GCGGITHASH),$(LAST_GCGGITHASH))
		@-$(MAKE) githash
endif
		@$(SHELL) -ec 'if test ! -e $(GCGGITHASHFILE) ; \
			then \
				echo "-> generating $(GCGGITHASHFILE)" ; \
				$(MAKE) githash ; \
			fi'
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_GCGGITHASH=$(GCGGITHASH)" >> $(LASTSETTINGS)
		@echo "LAST_LPS=$(LPS)" >> $(LASTSETTINGS)
		@echo "LAST_BLISS=$(BLISS)" >> $(LASTSETTINGS)
		@echo "LAST_CLIQUER=$(CLIQUER)" >> $(LASTSETTINGS)
		@echo "LAST_USRFLAGS=$(USRFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USROFLAGS=$(USROFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRDFLAGS=$(USRDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_OPENMP=$(OPENMP)" >> $(LASTSETTINGS)
		@echo "LAST_CPLEXSOLVER=$(CPLEXSOLVER)" >> $(LASTSETTINGS)
		@echo "LAST_STATISTICS=$(STATISTICS)" >> $(LASTSETTINGS)

.PHONY: $(SOFTLINKS)
$(SOFTLINKS): $(LIBDIR)/static $(LIBDIR)/shared $(LIBDIR)/include
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				echo ; \
		                echo "* GCG needs links to external libraries" ; \
		                echo "* Please insert the paths to the libaries below." ; \
		                echo "* The links will be installed in the 'lib' directory." ; \
		                echo ; \
		                echo -e $(LINKMSG) ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(GCGDIR) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					pwd;\
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

$(SCIPDIR): |$(LIBDIR)
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				echo ; \
		                echo "* GCG needs a link to SCIP" ; \
		                echo "* Please insert the paths to SCIP below." ; \
		                echo "* The link will be installed in the 'lib' directory." ; \
		                echo ; \
				echo "> Enter soft-link target file or directory for \"scip\" (e.g., scipoptsuite-4.0.0/scip-4.0.0):" ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(GCGDIR) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					pwd;\
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

.PHONY: help
help:
		@echo "Use the GCG makefile system."
		@echo
		@echo "  The main options for the GCG makefile system are as follows:"
		@echo
		@echo "  General options:"
		@echo "  - OPT={opt|dbg|prf}: Set solver mode (default: opt)."
		@echo "  - STATISTICS=<true|false>: Enable additional statistics. Required for pricing visualizations."
		@echo "  - MODE={readdec|none}: If set to readdec (default), GCG looks for given .dec files. "
		@echo "  - PARASCIP=<true|false>: Use SCIP's parallelization."
		@echo "  - OPENMP=<true|false>: Use GCG's parallelization. Will set PARASCIP to true."
		@echo
		@echo "  Additional Features and Modules:"
		@echo "  - READLINE=<true|false>: Enables READLINE, required for command line interaction (default: true)."
		@echo "  - CLIQUER=<true|false>: Enables CLIQUER (as a heuristic for stable set pricing problems)."
		@echo "  - HMETIS=<true|false>: Enables hMETIS (hypergraph partitioning, used in structure detection)."
		@echo "  - GSL=<true|false>: Enables the GNU Scientific Library (needed by a detector)"
		@echo "  - GAMS=<true|false>: To enable or disable (default) reading functionality in GAMS reader (needs GAMS)."
		@echo "  - GTEST=<true|false>: Enables Google Test."
		@echo "  - SYM=<none|bliss|sbliss|nauty|snauty>: To choose type of symmetry handling."
		@echo "  - ZIMPL=<true|false>: Enables ZIMPL, required to convert .zpl files to .lp/.mps files"
		@echo
		@echo "  More detailed options:"
		@echo "  - VALGRIND=<true|false>: Enable memory leak checking (and more) using valgrind."
		@echo "  - EXPRINT=<cppad|none>: Use CppAD as expressions interpreter (default) or no expressions interpreter."
		@echo "  - IPOPT=<true|false>: Turns support of IPOPT on or off (default)."
		@echo "  - LPSOPT=<dbg|opt>: Use debug or optimized (default) mode for LP-solver (SoPlex and Clp only)."
		@echo "  - NOBLKBUFMEM=<true|false>: Turn usage of internal memory functions off or on (default)."
		@echo "  - NOBLKMEM=<true|false>: Turn off block memory or on (default)."
		@echo "  - NOBUFMEM=<true|false>>: Turn off buffer memory or on (default)."
		@echo "  - ZIMPLOPT=<dbg|opt>: Use debug or optimized (default) mode for ZIMPL."
		@echo "  - ARCH"
		@echo "  - GMP"
		@echo "  - IPOPTOPT"
		@echo "  - LPSCHECK"
		@echo "  - PROJECT"
		@echo "  - SANITIZE"
		@echo "  - ZLIB"
		@echo
		@echo "  Options for make test:"
		@echo "  - TEST=file: Define a testset file (located in ./check/testset) to be used"
		@echo "  - SETTINGS=file: Define a settings file (located in ./settings) to be used"
		@echo "  - MASTERSETTINGS: Define a settings file for master problem(located in ./settings) to be used"
		@echo "  - MODE={readdec|none}: If set to readdec (default), GCG looks for given .dec files."
		@echo "  - STATISTICS=<true|false>: Enable additional statistics. Required for visualizations."
		@echo "  - OPT={opt|dbg|prf}: Set solver mode (default: opt)."
		@echo "  - MEM=b: Set memory limit."
		@echo "  - TIME=s: Set time limit in seconds."
		@echo "  - NODES=n: Set opened node limit for the branch and bound tree."
		@echo "  - DETECTIONSTATISTICS=<true|false>: Print extended detection statistics. Required for detection visualizations."
		@echo
		@echo "  Options for the Visualization Suite: (for targets, see last section)"
		@echo "  - VISU=<true|false>: Flag for make test. If true, generate data, then compose a testset report."
		@echo "  - DATADIR=folder: Directory including .out and .res files for report generation."
		@echo "  - VISUSETTINGS=file: Define a settings file (located in the root directory) for visualization scripts to use."
		@echo
		@echo "  Targets common for SCIP and GCG can be found in SCIP's make help."
		@echo "  Most important SCIP targets:"
		@echo "  - all (default): Build SCIP libaries and binary."
		@echo "  - links: Reconfigures the links in the \"lib\" directory."
		@echo "  - doc: Creates documentation in ./doc."
		@echo "  - check/test: Runs the check/test script current computer."
		@echo "  GCG specific targets:"
		@echo "  - deps: build all dependencies."
		@echo "  - testcluster: Runs the check/test script on the OR cluster."
		@echo "  - visu: create a PDF report from existing runtime data."


#---- EOF --------------------------------------------------------------------
