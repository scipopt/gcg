#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       *
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

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------
VERSION         :=	2.1.4
GCGGITHASH	=
SCIPDIR         =   lib/scip

#-----------------------------------------------------------------------------
# necessary information
#-----------------------------------------------------------------------------
LIBDIR          =	lib
DIRECTORIES     =	$(LIBDIR) $(LIBOBJDIR) $(LIBOBJSUBDIRS)
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
OPENMP      =   false
GSL         =   false
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings
LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(BLISS)

# overriding SCIP PARASCIP setting if compiled with OPENMP
ifeq ($(OPENMP),true)
override PARASCIP=true
endif

# overriding SCIP LPS setting if compiled with CPLEXSOLVER
ifeq ($(CPLEXSOLVER),true)
override LPS =  cpx
endif


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
-include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# BLISS
#-----------------------------------------------------------------------------

ifeq ($(BLISS),false)
FLAGS		+=	-DNBLISS
else
LDFLAGS		+= 	-lbliss
ifeq ($(COMP),gnu)
FLAGS		+=	-isystem$(LIBDIR)/blissinc
else
FLAGS		+=	-I$(LIBDIR)/blissinc
endif
SOFTLINKS	+=	$(LIBDIR)/blissinc
SOFTLINKS	+=	$(LIBDIR)/libbliss.$(STATICLIBEXT)
LINKMSG		+=	"bliss graph isomorphism framework (disable by compiling with \"make BLISS=false\"):\n"
LINKMSG		+=	" -> blissinc is the path to the bliss include files, e.g., \"bliss-0.72\"\n"
LINKMSG		+=	" -> \"libbliss.$(STATICLIBEXT)\" is the path to the bliss library, e.g., \"blissinc/libbliss.$(STATICLIBEXT)\"\n"
endif

#-----------------------------------------------------------------------------
# GSL
#-----------------------------------------------------------------------------
ifeq ($(GSL),true)
LDFLAGS                +=      -lgsl -lgslcblas -lm
FLAGS		           +=	   -DGSL
endif

#-----------------------------------------------------------------------------
# CPLEX pricing solver
#-----------------------------------------------------------------------------

ifeq ($(CPLEXSOLVER),true)
FLAGS		+=	-DCPLEXSOLVER -I$(SCIPDIR)/lib/include/cpxinc
else
FLAGS		+=	-DNCPLEXSOLVER
endif

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gcg

LIBOBJ		=	reader_blk.o \
			reader_dec.o \
			reader_ref.o \
			gcgplugins.o \
			relax_gcg.o \
			pricer_gcg.o \
			branch_orig.o \
			branch_ryanfoster.o \
			branch_generic.o \
			cons_origbranch.o \
			cons_masterbranch.o \
			cons_integralorig.o \
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
			heur_gcgsimplerounding.o \
			heur_gcgshifting.o \
			heur_gcgveclendiving.o \
			heur_gcgzirounding.o \
			heur_greedycolsel.o \
			heur_masterdiving.o \
			heur_mastercoefdiving.o \
			heur_masterfracdiving.o \
			heur_masterlinesdiving.o \
			heur_mastervecldiving.o \
			heur_origdiving.o \
			heur_relaxcolsel.o \
			heur_restmaster.o \
			heur_setcover.o \
			heur_xpcrossover.o \
			heur_xprins.o \
			branch_empty.o \
			branch_relpsprob.o \
			masterplugins.o \
			nodesel_master.o \
			sepa_master.o \
			sepa_basis.o \
			disp_gcg.o \
			disp_master.o \
			dialog_gcg.o \
			dialog_master.o \
			event_bestsol.o \
			event_mastersol.o \
			event_relaxsol.o \
			event_solvingstats.o \
			event_display.o \
			solver_mip.o \
			solver_knapsack.o \
			cons_decomp.o \
			decomp.o \
			dec_arrowheur.o \
			dec_stairheur.o \
			dec_connected.o \
			dec_cutpacking.o \
			dec_staircase.o \
			dec_random.o \
			dec_colors.o \
			gcggithash.o \
			reader_gp.o \
			scip_misc.o \
			misc.o \
			gcgheur.o \
			gcgvar.o \
			class_pricingtype.o \
			class_stabilization.o \
			graph/weights.o \
			graph/inst.o \
			graph/graph_tclique.o \
			stat.o \
			objdialog.o \
			dialog_graph.o \
			gcgpqueue.o \
			gcgcol.o \
			class_colpool.o

ifeq ($(BLISS),true)
LIBOBJ		+=	bliss_automorph.o \
			dec_isomorph.o \
			bliss.o
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
GCGLIBFILE	=	$(LIBDIR)/lib$(GCGLIB).$(LIBEXT)
GCGLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(GCGLIBOBJ))
GCGLIBSRC	=	$(filter $(wildcard $(SRCDIR)/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
GCGLIBSRC	+=	$(filter $(wildcard $(SRCDIR)/*/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard */*/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
GCGLIBDEP	=	$(SRCDIR)/depend.gcglib.$(OPT)
GCGLIBLINK	=	$(LIBDIR)/lib$(GCGLIBSHORTNAME).$(BASE).$(LIBEXT)
GCGLIBSHORTLINK = 	$(LIBDIR)/lib$(GCGLIBSHORTNAME).$(LIBEXT)

ALLSRC		=	$(MAINSRC) $(GCGLIBSRC)
LDFLAGS		+=	$(LINKCXX_L)$(LIBDIR)
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
		@$(MAKE) -C $(SCIPDIR) $^ libs

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
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(MASTERSETTINGS) $(notdir $(BINDIR)/$(GCGLIBNAME).$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(MODE);

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
		$(SHELL) -ec '$(DCC) $(subst isystem,I,$(FLAGS)) $(DFLAGS) $(GCGLIBSRC) \
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

$(GCGLIBFILE):	$(GCGLIBOBJFILES)
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
ifneq ($(LAST_BLISS),$(BLISS))
		@-touch $(SRCDIR)/dec_isomorph.cpp
		@-touch $(SRCDIR)/relax_gcg.c
		@-touch $(SRCDIR)/gcgplugins.c
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
		@echo "LAST_USRFLAGS=$(USRFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USROFLAGS=$(USROFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRDFLAGS=$(USRDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_OPENMP=$(OPENMP)" >> $(LASTSETTINGS)
		@echo "LAST_CPLEXSOLVER=$(CPLEXSOLVER)" >> $(LASTSETTINGS)

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
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
