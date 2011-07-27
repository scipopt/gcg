#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Colum Generation                                 *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

#@file    Makefile
#@brief   Makefile for generic column generation code using SCIP as a callable library
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Marc Pfetsch
#@author  Gerald Gamrath


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       lib/scip


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project




#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gcg
MAINOBJ		=	reader_blk.o \
         		reader_ref.o \
			interface.o \
			gcgplugins.o \
			relax_gcg.o \
			pricer_gcg.o \
			branch_orig.o \
			branch_ryanfoster.o \
			cons_origbranch.o \
			cons_masterbranch.o \
			cons_integralOrig.o \
			heur_feasrestmaster.o \
			heur_gcgcoefdiving.o \
			heur_gcgfracdiving.o \
			heur_gcgrens.o \
			heur_gcgrounding.o \
			heur_gcgsimplerounding.o \
			heur_gcgshifting.o \
			heur_gcgzirounding.o \
			heur_greedycolsel.o \
			heur_relaxcolsel.o \
			heur_restmaster.o \
			branch_master.o \
			branch_relpsprob.o \
			masterplugins.o \
			pricingplugins.o \
			nodesel_master.o \
			sepa_master.o \
			disp_gcg.o \
			disp_master.o \
			dialog_gcg.o \
			solver_mip.o \
			solver_knapsack.o \
			cons_decomp.o \
			dec_arrowheur.o \
			dec_stairheur.o \
			dec_borderheur.o \
			reader_gp.o \
			scip_misc.o \
			main.o

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: test
test:		
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(BINDIR)/gcg.$(BASE).$(LPS) $(SETTINGS) $(notdir $(BINDIR)/gcg.$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK);

.PHONY: eval
eval:	
		cd check; \
		$(SHELL) ./eval.sh $(TEST) $(BINDIR)/gcg.$(BASE).$(LPS) $(SETTINGS) $(notdir $(BINDIR)/gcg.$(BASE).$(LPS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK);

.PHONY: clean
clean:
ifneq ($(OBJDIR),)
		-rm -rf $(OBJDIR)/*
endif
		-rm -f $(MAINFILE)

.PHONY: tags
tags:
		cd src/; rm -f TAGS; etags *.c *.h ../$(SCIPDIR)/src/scip/*.c ../$(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
