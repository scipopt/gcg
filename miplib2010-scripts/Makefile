#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.19 2011/06/04 20:34:24 bzfheinz Exp $

VERSION		=	1.0.2
TIME     	=  	3600
TEST		=	benchmark
SOLVER		=	scip
HARDMEM		=	8192
THREADS		= 	0

SHELL		= 	bash
DOXY		=	doxygen

CHECKERDIR      =       checker

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.PHONY: help
help:
		@echo "See README for details about the MIPLIB2010 test environment"
		@echo
		@echo "VERSION:      $(VERSION)"
		@echo
		@echo "TARGETS:"
		@echo "** checker -> compiles the solution checker" 
		@echo "** clean   -> cleans the solution checker" 
		@echo "** cmpres  -> generates solver comparison file"
		@echo "** doc     -> generates doxygen documentation"
		@echo "** eval    -> evaluate test run" 
		@echo "** test    -> start automatic test runs" 
		@echo 
		@echo "PARAMETERS:"
		@echo "** HARDMEM -> maximum memory to use MB [8192]"
		@echo "** SOLVER  -> solver [scip]"
		@echo "** THREADS -> number of threads (0: automatic) [0]"
		@echo "** TIME    -> time limit per instance in seconds [3600]"
		@echo "** TEST    -> tesr set [benchmark]"

.PHONY: checker
checker:
		@$(MAKE) -C $(CHECKERDIR) $^

.PHONY:		clean
clean: 		
		@$(MAKE) -C $(CHECKERDIR) clean

.PHONY: doc
doc: 		
		cd doc; $(DOXY) miplib.dxy;

.PHONY: test
test:
		@echo "run test: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST) TIME=$(TIME) HARDMEN=$(HARDMEM) THREADS=$(THREADS)"
		@$(SHELL) ./scripts/run.sh $(SHELL) $(VERSION) $(SOLVER) $(TEST) $(TIME) $(HARDMEM) $(THREADS);

.PHONY: eval
eval:
		@echo "evaluate test: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST)"
		@$(SHELL) ./scripts/evalrun.sh results/$(TEST).$(SOLVER).out;

.PHONY: cmpres
cmpres:
		@echo "compare result tables: VERSION=$(VERSION) SOLVER=$(SOLVER) TEST=$(TEST)"
		@$(SHELL) ./scripts/allcmpres.sh results/$(TEST).$(SOLVER).res;


# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
