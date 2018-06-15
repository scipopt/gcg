#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       *
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
VERSION         :=	3.0.0
GCGGITHASH	=
SCIPDIR         =   lib/scip

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

VALGRIND	=	false
MODE		=	readdec
STATISTICS  =  false
PROJECT		=	none
GTEST		=	false
PARASCIP	= 	true
BLISS      	=   false
CLIQUER     =   false
HMETIS      =   true
OPENMP      =   false
GSL         =   false
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings
LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(BLISS).$(CLIQUER)

# overriding SCIP PARASCIP setting if compiled with OPENMP
ifeq ($(OPENMP),true)
override PARASCIP=true
endif

# overriding SCIP LPS setting if compiled with CPLEXSOLVER
ifeq ($(CPLEXSOLVER),true)
override LPS =  cpx
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
			branch_empty.o \
			branch_generic.o \
			benders_gcg.o \
			branch_orig.o \
			branch_relpsprob.o \
			branch_ryanfoster.o \
			class_consclassifier.o \
			class_indexclassifier.o \
			class_miscvisualization.o \
			class_pricingcontroller.o \
			class_pricingtype.o \
			class_seeed.o \
			class_seeedpool.o \
			class_stabilization.o \
			class_varclassifier.o \
			colpool.o \
			cons_decomp.o \
			cons_integralorig.o \
			cons_masterbranch.o \
			cons_origbranch.o \
			dec_colors.o \
			dec_compgreedily.o \
			dec_connected.o \
			dec_connected_noNewLinkingVars.o \
			dec_connectedbase.o \
			dec_consclass.o \
			dec_consname.o \
			dec_constype.o \
			dec_cutpacking.o \
			dec_dbscan.o \
			dec_densemasterconss.o \
			dec_generalmastersetcover.o \
			dec_generalmastersetpack.o \
			dec_generalmastersetpart.o \
			dec_hcgpartition.o \
			dec_hrcgpartition.o \
			dec_hrgpartition.o \
			dec_mastersetcover.o \
			dec_mastersetpack.o \
			dec_mastersetpart.o \
			dec_mcl.o \
			dec_mst.o \
			dec_postprocess.o \
			dec_random.o \
			dec_staircase.o \
			dec_staircase_lsp.o \
			dec_stairheur.o \
			dec_varclass.o \
			decomp.o \
			dialog_gcg.o \
			dialog_graph.o \
			dialog_master.o \
			disp_gcg.o \
			disp_master.o \
			event_bestsol.o \
			event_display.o \
			event_mastersol.o \
			event_relaxsol.o \
			event_solvingstats.o \
			gcgcol.o \
			gcggithash.o \
			gcgheur.o \
			gcgsepa.o \
			gcgplugins.o \
			gcgpqueue.o \
			gcgsort.o \
			gcgvar.o \
			graph/graph_gcg.o \
			graph/graph_tclique.o \
			graph/inst.o \
			graph/weights.o \
			heur_gcgcoefdiving.o \
			heur_gcgdins.o \
			heur_gcgfeaspump.o \
			heur_gcgfracdiving.o \
			heur_gcgguideddiving.o \
			heur_gcglinesdiving.o \
			heur_gcgpscostdiving.o \
			heur_gcgrens.o \
			heur_gcgrins.o \
			heur_gcgrounding.o \
			heur_gcgshifting.o \
			heur_gcgsimplerounding.o \
			heur_gcgveclendiving.o \
			heur_gcgzirounding.o \
			heur_greedycolsel.o \
			heur_mastercoefdiving.o \
			heur_masterdiving.o \
			heur_masterfracdiving.o \
			heur_masterlinesdiving.o \
			heur_mastervecldiving.o \
			heur_origdiving.o \
			heur_relaxcolsel.o \
			heur_restmaster.o \
			heur_setcover.o \
			heur_xpcrossover.o \
			heur_xprins.o \
			masterplugins.o \
			bendersplugins.o \
			misc.o \
			nodesel_master.o \
			objdialog.o \
			params_visu.o \
			presol_roundbound.o \
			pricer_gcg.o \
			pricestore_gcg.o \
			pricingjob.o \
			pricingprob.o \
			reader_blk.o \
			reader_cls.o \
			reader_dec.o \
			reader_gp.o \
			reader_ref.o \
			reader_tex.o \
			relax_gcg.o \
			scip_misc.o \
			sepa_basis.o \
			sepa_master.o \
			solver.o \
			solver_knapsack.o \
			solver_mip.o \
			stat.o \

ifeq ($(BLISS),true)
LIBOBJ		+=	bliss_automorph.o \
			dec_isomorph.o \
			bliss.o
endif

ifeq ($(CLIQUER),true)
LIBOBJ		+=	solver_cliquer.o
endif

ifeq ($(CPLEXSOLVER),true)
LIBOBJ		+=	solver_cplex.o
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
OBJSUBDIRS	= 	graph
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
GCGLIBSRC	+=	$(filter $(wildcard $(SRCDIR)/*/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard */*/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
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
		cd doc; $(DOXY) $(MAINNAME).dxy;

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

.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(MASTERSETTINGS) $(notdir $(BINDIR)/$(GCGLIBNAME).$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(MODE) $(SETCUTOFF) $(STATISTICS);

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


$(MAINFILE):	$(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(TPILIBFILE) $(MAINOBJFILES) $(GCGLIBFILE) | $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_l)$(GCGLIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_L)$(SCIPDIR)/lib/static $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(TPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.c | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@

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
		@-touch $(SRCDIR)/pricer_gcg.h
		@-touch $(SRCDIR)/pricer_gcg.cpp
		@-touch $(SRCDIR)/stat.c
		@-touch $(SRCDIR)/event_bestsol.h
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


#---- EOF --------------------------------------------------------------------
