#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
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
VERSION         :=	1.0.0.1
SCIPDIR         =       lib/scip
#-----------------------------------------------------------------------------
# necessary information
#-----------------------------------------------------------------------------


LIBDIR          =       lib
DIRECTORIES     =       $(LIBDIR)
MAKESOFTLINKS	=	true

SHELL		= 	bash
READ		=	read -e
LN_s		= 	ln -s
GCGDIR		=	$(realpath .)

VALGRIND        =       false
DECMODE		=	readdec

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
-include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gcg
TESTNAME	=	gcg_test
LIBOBJ		=	reader_blk.o \
			reader_dec.o \
			reader_ref.o \
			gcgplugins.o \
			relax_gcg.o \
			pricer_gcg.o \
			branch_orig.o \
			branch_ryanfoster.o \
			cons_origbranch.o \
			cons_masterbranch.o \
			cons_integralorig.o \
			heur_gcgcoefdiving.o \
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
			heur_relaxcolsel.o \
			heur_restmaster.o \
			heur_xpcrossover.o \
			heur_xprins.o \
			branch_master.o \
			branch_relpsprob.o \
			masterplugins.o \
			pricingplugins.o \
			nodesel_master.o \
			sepa_master.o \
			disp_gcg.o \
			disp_master.o \
			dialog_gcg.o \
			dialog_master.o \
			event_solvingstats.o \
			solver_mip.o \
			solver_knapsack.o \
			cons_decomp.o \
			decomp.o \
			dec_arrowheur.o \
			dec_borderheur.o \
			dec_stairheur.o \
			dec_connected.o \
			dec_cutpacking.o \
			dec_staircase.o \
			dec_random.o \
			gcggithash.o \
			reader_gp.o \
			scip_misc.o \
			misc.o \
			gcgvar.o \
			class_pricingtype.o \
			stat.o

MAINOBJ		=	${LIBOBJ} \
			main.o

TESTOBJ		=	tests/test.o

MAINSRC		=	$(filter $(wildcard $(SRCDIR)/*.c),$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*.cpp),$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp)))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

TESTSRC		=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.c)) $(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.cpp))
TESTDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

TEST		=	$(TESTNAME).$(BASE).$(LPS)$(EXEEXTENSION)
TESTFILE	=	$(BINDIR)/$(TEST)
TESTSHORTLINK	=	$(BINDIR)/$(TESTNAME)
TESTOBJFILES	=	$(addprefix $(OBJDIR)/,$(TESTOBJ))
TESTOBJDIR      =       $(OBJDIR)/tests

CFLAGS          +=       -Ilib/gtest/
CXXFLAGS        +=       -Ilib/gtest/
LDFLAGS         +=       -Llib/ -lgtest


SOFTLINKS	+=	$(LIBDIR)/scip
LPIINSTMSG	=	"  -> \"scip\" is the path to the SCIP directory, e.g., \"scipoptsuite-3.0.0/scip-3.0.0/\""
LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.scip

# GCG Library
LIBOBJDIR	=	$(OBJDIR)/lib
LIBOBJSUBDIRS	=       

GCGLIBSHORTNAME =	gcg
GCGLIBNAME	=	$(GCGLIBSHORTNAME)-$(VERSION)

GCGLIBOBJ	=	${LIBOBJ}
GCGLIB		=	$(GCGLIBNAME).$(BASE)
GCGLIBFILE	=	$(LIBDIR)/lib$(GCGLIB).$(LIBEXT)
GCGLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(GCGLIBOBJ))
GCGLIBSRC	=	$(filter $(wildcard $(SRCDIR)/*.c),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.c))) $(filter $(wildcard $(SRCDIR)/*.cpp),$(addprefix $(SRCDIR)/,$(GCGLIBOBJ:.o=.cpp)))
GCGLIBDEP	=	$(SRCDIR)/depend.gcglib.$(OPT)
GCGLIBLINK	=	$(LIBDIR)/lib$(GCGLIBSHORTNAME).$(BASE).$(LIBEXT)
GCGLIBSHORTLINK = 	$(LIBDIR)/lib$(GCGLIBSHORTNAME).$(LIBEXT)



#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------


ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK) ${GCGLIBFILE} ${GCGLIB} ${GCGLIBSHORTLINK} ${TESTSHORTLINK} ${LIBOBJFILES} ${TESTFILE} ${TESTMAIN}
endif

ifeq ($(OPENMP),true)
CFLAGS+="-fopenmp"
LDFLAGS+="-fopenmp"
CXXFLAGS+="-fopenmp"
endif

.PHONY: all
all:       githash $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

$(SCIPDIR)/make/make.project: $(LINKSMARKERFILE);

.PHONY: libs
libs:		$(GCGLIBFILE) $(GCGLIBSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
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


.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy; cp tabs.css html/

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

# include target to detect the current git hash
-include make/local/make.detectgithash
-include make/local/make.targets

# this empty target is needed for the SCIP release versions
githash::   # do not remove the double-colon

.PHONY: test
test:
		cd check; \
		echo $(SHELL) ./check.sh $(TEST) $(BINDIR)/gcg.$(BASE).$(LPS) $(SETTINGS) $(notdir $(BINDIR)/gcg.$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(DECMODE); \
		$(SHELL) ./check.sh $(TEST) $(BINDIR)/gcg.$(BASE).$(LPS) $(SETTINGS) $(notdir $(BINDIR)/gcg.$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(DECMODE);

.PHONY: tests
tests: 		libs $(TESTFILE) $(TESTSHORTLINK)

$(TESTSHORTLINK):	$(TESTFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(TESTFILE)) $(notdir $@)

$(TESTFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(GCGLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(TESTOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(TESTOBJFILES)  $(LINKCXX_l)$(GCGLIBSHORTNAME)$(LINKLIBSUFFIX) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX)  \
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS)  -L$(LIBDIR)/gtest/lib/ -L$(LIBDIR) -lgtest $(LINKCXX_o)$@


.PHONY: eval
eval:
		cd check; \
		$(SHELL) ./eval.sh $(TEST) $(BINDIR)/gcg.$(BASE).$(LPS) $(SETTINGS) $(notdir $(BINDIR)/gcg.$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND);

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		-(cd ./$(LIBOBJDIR) && rm -f *.o)
		-rmdir $(LIBOBJDIR)
		-(cd ./$(TESTOBJDIR) && rm -f *.o)
		-rmdir $(TESTOBJDIR)
		-(cd ./$(OBJDIR) && rm -f *.o)
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE) $(MAINSHORTLINK)

.PHONY: tags
tags:
		cd src/; rm -f TAGS; etags *.c *.h ../$(SCIPDIR)/src/scip/*.c ../$(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR) gcgdepend
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'
-include	$(MAINDEP)

.PHONY: gcgdepend
gcgdepend:
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(GCGLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(GCGLIBDEP)'

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@


$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

$(GCGLIBFILE):	$(LIBOBJDIR) $(LIBDIR) $(LIBOBJSUBDIRS)  $(GCGLIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(GCGLIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(LIBOBJDIR):	$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(TESTOBJDIR):	$(OBJDIR)
		@-mkdir -p $(TESTOBJDIR)

$(GCGLIBSHORTLINK):	$(GCGLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(GCGLIBFILE)) $(notdir $@)


$(LINKSMARKERFILE): links
#		@$(MAKE) links

.PHONY: links
links:		$(LIBDIR) $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

$(DIRECTORIES):
		@echo
		@echo "- creating directory \"$@\""
		@-mkdir -p $@

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				echo ; \
		                echo "* GCG needs a softlink to SCIP" ; \
		                echo "* Please insert the paths to scip below." ; \
		                echo "* The link will be installed in the 'lib' directory." ; \
		                echo "* For more information and if you experience problems see the INSTALL file." ; \
		                echo ; \
		                echo -e $(LPIINSTMSG) ; \
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


#---- EOF --------------------------------------------------------------------
