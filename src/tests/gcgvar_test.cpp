/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2020 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gcgvar_test.cpp
 * @brief  Description
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "test.h"
#include "pub_gcgvar.h"
#include "struct_vardata.h"
#include "scip/struct_var.h"
#include "relax_gcg.h"
#include "gcg.h"

#define ORIGVAR(ovar, ovardata) SCIP_VAR ovar;  SCIP_VARDATA ovardata; ovar.vardata = &ovardata;   ovardata.vartype = GCG_VARTYPE_ORIGINAL
#define PRICINGVAR(pvar, pvardata) SCIP_VAR pvar;  SCIP_VARDATA pvardata; pvar.vardata =&pvardata;   pvardata.vartype = GCG_VARTYPE_PRICING
#define MASTERVAR(mvar, mvardata) SCIP_VAR mvar;  SCIP_VARDATA mvardata; mvar.vardata = &mvardata;   mvardata.vartype = GCG_VARTYPE_MASTER
#define LINKINGVAR(lvar, lvardata) SCIP_VAR lvar;  SCIP_VARDATA lvardata; lvar.vardata = &lvardata; lvardata.vartype = GCG_VARTYPE_ORIGINAL; lvardata.blocknr = -2

class GcgVarTest : public ::testing::Test {

   virtual void SetUp() {
      decomp = NULL;
      SCIP_CALL_ABORT( SCIPcreate(&scip) );
   }

   virtual void TearDown() {
      SCIP_CALL_ABORT( SCIPfree(&scip) );
   }
protected:
   DEC_DECOMP* decomp;

public:
   static SCIP* scip;
};

SCIP* GcgVarTest::scip = NULL;

TEST_F(GcgVarTest, PricingVarIsPricingVar) {
   PRICINGVAR(var, vardata);
   ASSERT_EQ(TRUE, GCGvarIsPricing(&var));
}

TEST_F(GcgVarTest, MasterVarIsNotPricingVar) {
   MASTERVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsPricing(&var));
}

TEST_F(GcgVarTest, OriginalVarIsNotPricingVar) {
   ORIGVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsPricing(&var));
}

TEST_F(GcgVarTest, PricingVarIsNotMasterVar) {
   PRICINGVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsMaster(&var));
}

TEST_F(GcgVarTest, MasterVarIsMasterVar) {
   MASTERVAR(var, vardata);
   ASSERT_EQ(TRUE, GCGvarIsMaster(&var));
}

TEST_F(GcgVarTest, OriginalVarIsNotMasterVar) {
   ORIGVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsPricing(&var));
}

TEST_F(GcgVarTest, PricingVarIsNotOriginalVar) {
   PRICINGVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsOriginal(&var));
}

TEST_F(GcgVarTest, MasterVarIsNotOriginalVar) {
   MASTERVAR(var, vardata);
   ASSERT_EQ(FALSE, GCGvarIsOriginal(&var));
}

TEST_F(GcgVarTest, OriginalVarIsOriginalVar) {
   ORIGVAR(var, vardata);
   ASSERT_EQ(TRUE, GCGvarIsOriginal(&var));
}

TEST_F(GcgVarTest, LinkingVarIsLinkingVar) {
   LINKINGVAR(var, vardata);
   ASSERT_EQ(TRUE, GCGoriginalVarIsLinking(&var));
}

TEST_F(GcgVarTest, BlockVarIsNotLinkingVar) {
   SCIP_VAR var;
   SCIP_VARDATA vardata;
   var.vardata = &vardata;
   vardata.blocknr = 1;
   vardata.vartype = GCG_VARTYPE_ORIGINAL;

   ASSERT_EQ(FALSE, GCGoriginalVarIsLinking(&var));
}

TEST_F(GcgVarTest, MasterVarIsNotLinkingVar) {
   SCIP_VAR var;
   SCIP_VARDATA vardata;
   var.vardata = &vardata;
   vardata.blocknr = -1;
   vardata.vartype = GCG_VARTYPE_ORIGINAL;

   ASSERT_EQ(FALSE, GCGoriginalVarIsLinking(&var));
}

TEST_F(GcgVarTest, OriginalVarGetPricingVar) {
   ORIGVAR(var, vardata);
   PRICINGVAR(pricingvar, pricingvardata);

   vardata.blocknr = 0;
   vardata.data.origvardata.pricingvar = &pricingvar;
   vardata.data.origvardata.linkingvardata = NULL;
   ASSERT_EQ(&pricingvar, GCGoriginalVarGetPricingVar(&var));
}

TEST_F(GcgVarTest, OriginalVarSetPricingVar) {
   ORIGVAR(var, vardata);
   PRICINGVAR(pricingvar, pricingvardata);

   vardata.blocknr = 0;

   vardata.data.origvardata.linkingvardata = NULL;
   GCGoriginalVarSetPricingVar(&var, &pricingvar);
   ASSERT_EQ(&pricingvar, vardata.data.origvardata.pricingvar);
}


TEST_F(GcgVarTest, LinkingVarGetPricingVars) {
   LINKINGVAR(var, vardata);
   GCG_LINKINGVARDATA linkvardata;

   vardata.data.origvardata.linkingvardata = &linkvardata;
   vardata.data.origvardata.linkingvardata->pricingvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ((SCIP_VAR**) 0xDEADBEEF, GCGlinkingVarGetPricingVars(&var));
}


TEST_F(GcgVarTest, LinkingVarSetPricingVar) {
   LINKINGVAR(lvar, lvardata);
   PRICINGVAR(pvar, pvardata);
   GCG_LINKINGVARDATA linkvardata;
   SCIP_VAR* vars[4] = {NULL, NULL, NULL, NULL};
   lvardata.data.origvardata.linkingvardata = &linkvardata;
   linkvardata.pricingvars = vars;
   GCGlinkingVarSetPricingVar(&lvar, 2, &pvar);
   ASSERT_EQ(NULL, vars[0]);
   ASSERT_EQ(NULL, vars[1]);
   ASSERT_EQ(&pvar, vars[2]);
   ASSERT_EQ(NULL, vars[3]);
}

TEST_F(GcgVarTest, LinkingVarGetBlocksArrayLargeEnough) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkvardata;
   SCIP_VAR* vars[4] = {NULL, (SCIP_VAR*)0xDEADBEEF, NULL, (SCIP_VAR*)0xDEADBEEF};
   int blocks[4] = {-1,-1,-1,-1};
   lvardata.data.origvardata.linkingvardata = &linkvardata;
   linkvardata.pricingvars = vars;
   linkvardata.nblocks = 2;

   SCIP_CALL_EXPECT(GCGlinkingVarGetBlocks(&lvar, 4, blocks));

   ASSERT_EQ(1, blocks[0]);
   ASSERT_EQ(3, blocks[1]);
   ASSERT_EQ(-1, blocks[2]);
   ASSERT_EQ(-1, blocks[3]);
}

TEST_F(GcgVarTest, LinkingVarGetBlocksArrayTooSmall) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkvardata;
   lvardata.data.origvardata.linkingvardata = &linkvardata;
   linkvardata.nblocks = 2;
   ASSERT_EQ(SCIP_INVALIDDATA, GCGlinkingVarGetBlocks(&lvar, 1, (int*) 0xDEADBEEF));
}

TEST_F(GcgVarTest, LinkingVarGetNBlocks) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkvardata;
   lvardata.data.origvardata.linkingvardata = &linkvardata;
   linkvardata.nblocks = 4;
   ASSERT_EQ(4, GCGlinkingVarGetNBlocks(&lvar));
}

TEST_F(GcgVarTest, PricingVarGetOriginalVar) {
   PRICINGVAR(var, vardata);
   ORIGVAR(ovar, ovardata);
   SCIP_VAR* vars[1] = {&ovar};
   vardata.blocknr = 0;
   vardata.data.pricingvardata.origvars = vars;
   vardata.data.pricingvardata.norigvars = 1;
   ASSERT_EQ(&ovar, GCGpricingVarGetOriginalVar(&var));
}

TEST_F(GcgVarTest, PricingVarGetOrigvars) {
   PRICINGVAR(var, vardata);
   ORIGVAR(ovar, ovardata);
   SCIP_VAR* vars[1] = {&ovar};
   vardata.blocknr = 0;
   vardata.data.pricingvardata.origvars = vars;
   vardata.data.pricingvardata.norigvars = 1;
   ASSERT_EQ(vars, GCGpricingVarGetOrigvars(&var));
}

TEST_F(GcgVarTest, PricingVarGetNOriginalVars) {
   PRICINGVAR(var, vardata);

   SCIP_VAR* vars[1]= {NULL};
   vardata.blocknr = 0;
   vardata.data.pricingvardata.origvars = vars;
   vardata.data.pricingvardata.norigvars = 1;
   ASSERT_EQ(1, GCGpricingVarGetNOrigvars(&var));
}

TEST_F(GcgVarTest, PricingVarAddOriginalVarWhenNonempty) {
   PRICINGVAR(var, vardata);
   ORIGVAR(ovar, ovardata);

   SCIP_VAR** vars;
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &vars, 1));
   vars[0]= &ovar;

   vardata.blocknr = 0;
   vardata.data.pricingvardata.origvars = vars;
   vardata.data.pricingvardata.norigvars = 1;
   vardata.data.pricingvardata.maxorigvars = 1;

   SCIP_CALL_EXPECT(GCGpricingVarAddOrigVar(scip, &var, &ovar));
   ASSERT_EQ(2, GCGpricingVarGetNOrigvars(&var));
   SCIPfreeBlockMemoryArray(scip, &vardata.data.pricingvardata.origvars, vardata.data.pricingvardata.maxorigvars);
}

TEST_F(GcgVarTest, PricingVarAddOriginalVarWhenEmpty) {
   PRICINGVAR(var, vardata);
   ORIGVAR(ovar, ovardata);

   SCIP_VAR** vars;
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &vars, 1));
   vars[0]= (SCIP_VAR*)0xDEADBEEF;

   vardata.blocknr = 0;
   vardata.data.pricingvardata.origvars = vars;
   vardata.data.pricingvardata.norigvars = 0;
   vardata.data.pricingvardata.maxorigvars = 1;

   SCIP_CALL_EXPECT(GCGpricingVarAddOrigVar(scip, &var, &ovar));
   ASSERT_EQ(1, GCGpricingVarGetNOrigvars(&var));
   SCIPfreeBlockMemoryArray(scip, &vardata.data.pricingvardata.origvars, vardata.data.pricingvardata.maxorigvars);
}

TEST_F(GcgVarTest, OriginalVarGetNMastervars) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.nmastervars = 0xDEAD;
   ASSERT_EQ(0xDEAD, GCGoriginalVarGetNMastervars(&ovar));
}

TEST_F(GcgVarTest, OriginalVarGetMastervars) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.mastervars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ((SCIP_VAR**)0xDEADBEEF, GCGoriginalVarGetMastervars(&ovar));
}

TEST_F(GcgVarTest, OriginalVarGetMastervals) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.mastervals = (SCIP_Real*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_Real*)0xDEADBEEF, GCGoriginalVarGetMastervals(&ovar));
}

TEST_F(GcgVarTest, OriginalVarGetCoefs) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.coefs = (SCIP_Real*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_Real*)0xDEADBEEF, GCGoriginalVarGetCoefs(&ovar));
}

TEST_F(GcgVarTest, OriginalVarGetNCoefs) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.coefs = (SCIP_Real*) 0xDEADBEEF;
   ovardata.data.origvardata.ncoefs = 0xDEAD;
   ASSERT_EQ(0xDEAD, GCGoriginalVarGetNCoefs(&ovar));
}

TEST_F(GcgVarTest, OriginalVarSetNCoefs) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.ncoefs = 0;
   ovardata.data.origvardata.coefs = (SCIP_Real*) 0xDEADBEEF;
   GCGoriginalVarSetNCoefs(&ovar, 0xDEAD);
   ASSERT_EQ(0xDEAD, GCGoriginalVarGetNCoefs(&ovar));
}

TEST_F(GcgVarTest, OriginalVarAddCoefsWhenEmpty) {
   ORIGVAR(ovar, ovardata);
   SCIP_CONS* cons = (SCIP_CONS*) 0xDEADBEEF;
   ovardata.data.origvardata.ncoefs = 0;
   ovardata.data.origvardata.coefs = NULL;
   ovardata.blocknr = 0;
   ovardata.data.origvardata.masterconss = NULL;

   SCIP_CALL_EXPECT(GCGoriginalVarAddCoef(scip, &ovar, 1.0, cons));
   ASSERT_EQ(1, ovardata.data.origvardata.ncoefs);
   ASSERT_EQ(1.0, ovardata.data.origvardata.coefs[0]);
   ASSERT_EQ(cons, ovardata.data.origvardata.masterconss[0]);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.coefs, ovardata.data.origvardata.ncoefs);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.masterconss, ovardata.data.origvardata.ncoefs);
}

TEST_F(GcgVarTest, OriginalVarAddCoefsWhenNonemty) {
   ORIGVAR(ovar, ovardata);
   SCIP_CONS* cons = (SCIP_CONS*) 0xDEADBEEF;
   ovardata.data.origvardata.ncoefs = 1;
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.coefs, 1));
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.masterconss, 1));
   ovardata.data.origvardata.coefs[0] = 1.0;
   ovardata.data.origvardata.masterconss[0] = (SCIP_CONS*) 0xDEADCAFF;

   SCIP_CALL_EXPECT(GCGoriginalVarAddCoef(scip, &ovar, 2.0, cons));
   ASSERT_EQ(2, ovardata.data.origvardata.ncoefs);
   ASSERT_EQ(1.0, ovardata.data.origvardata.coefs[0]);
   ASSERT_EQ((SCIP_CONS*) 0xDEADCAFF, ovardata.data.origvardata.masterconss[0]);
   ASSERT_EQ(2.0, ovardata.data.origvardata.coefs[1]);
   ASSERT_EQ(cons, ovardata.data.origvardata.masterconss[1]);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.coefs, ovardata.data.origvardata.ncoefs);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.masterconss, ovardata.data.origvardata.ncoefs);
}

TEST_F(GcgVarTest, OriginalVarGetMasterconss) {
   ORIGVAR(ovar, ovardata);
   ovardata.data.origvardata.masterconss = (SCIP_CONS**) 0xDEADBEEF;
   ASSERT_EQ((SCIP_CONS**)0xDEADBEEF, GCGoriginalVarGetMasterconss(&ovar));
}

TEST_F(GcgVarTest, OriginalVarAddFirstBlock) {
   ORIGVAR(ovar, ovardata);
   ovardata.blocknr = 0;

   SCIP_CALL_EXPECT(GCGoriginalVarAddBlock(scip, &ovar, 2, 4, DEC_DECMODE_DANTZIGWOLFE));
   ASSERT_NE((GCG_LINKINGVARDATA*) NULL, ovardata.data.origvardata.linkingvardata);

   ASSERT_NE((SCIP_VAR**) NULL, ovardata.data.origvardata.linkingvardata->pricingvars);
   ASSERT_NE((SCIP_CONS**) NULL, ovardata.data.origvardata.linkingvardata->linkconss);


   ASSERT_EQ(&ovar, ovardata.data.origvardata.linkingvardata->pricingvars[0]);
   ASSERT_EQ(&ovar, ovardata.data.origvardata.linkingvardata->pricingvars[2]);
   ASSERT_EQ(2, ovardata.data.origvardata.linkingvardata->nblocks);
   ASSERT_EQ(TRUE, GCGoriginalVarIsLinking(&ovar));

   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.linkingvardata->pricingvars, 4);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.linkingvardata->linkconss, 4);
   SCIPfreeBlockMemory(scip, &ovardata.data.origvardata.linkingvardata);

}

TEST_F(GcgVarTest, OriginalVarAddSecondBlock) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkingvardata;
   SCIP_VAR* pricingvars[3] = {&lvar, NULL, &lvar};
   SCIP_CONS* conss[3] = {(SCIP_CONS*) 0xDEADBEEF, NULL, (SCIP_CONS*) 0xDEADCAFF};
   linkingvardata.pricingvars = pricingvars;
   linkingvardata.nblocks = 2;
   linkingvardata.linkconss = conss;
   lvardata.data.origvardata.linkingvardata = &linkingvardata;
   SCIP_CALL_EXPECT(GCGoriginalVarAddBlock(scip, &lvar, 1, 3, DEC_DECMODE_DANTZIGWOLFE));
   ASSERT_EQ(&lvar, linkingvardata.pricingvars[0]);
   ASSERT_EQ(&lvar, linkingvardata.pricingvars[1]);
   ASSERT_EQ(&lvar, linkingvardata.pricingvars[2]);
   ASSERT_EQ(3, linkingvardata.nblocks);
}


TEST_F(GcgVarTest, VarGetBlock) {
   ORIGVAR(ovar, ovardata);
   ovardata.blocknr = 2;
   ASSERT_EQ(2, GCGvarGetBlock(&ovar));
}

TEST_F(GcgVarTest, VarSetBlock) {
   ORIGVAR(ovar, ovardata);
   ovardata.blocknr = 0;
   GCGvarSetBlock(&ovar, 2);
   ASSERT_EQ(2, ovardata.blocknr);
}

TEST_F(GcgVarTest, LinkingVarGetLinkingConss) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkingvardata;
   SCIP_CONS* conss[1] = {NULL};
   lvardata.data.origvardata.linkingvardata = &linkingvardata;
   linkingvardata.linkconss = conss;

   ASSERT_EQ(conss, GCGlinkingVarGetLinkingConss(&lvar));
}

TEST_F(GcgVarTest, LinkingVarSetLinkingConss) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkingvardata;
   SCIP_CONS* conss[3] = {(SCIP_CONS*) 0xDEADBEEF, NULL, (SCIP_CONS*) 0xDEADCAFF};
   SCIP_CONS* cons = (SCIP_CONS*) 0xDEADDEAD;

   lvardata.data.origvardata.linkingvardata = &linkingvardata;
   linkingvardata.linkconss = conss;

   GCGlinkingVarSetLinkingCons(&lvar, cons, 1);

   ASSERT_EQ((SCIP_CONS*) 0xDEADBEEF, linkingvardata.linkconss[0]);
   ASSERT_EQ((SCIP_CONS*) 0xDEADDEAD, linkingvardata.linkconss[1]);
   ASSERT_EQ((SCIP_CONS*) 0xDEADCAFF, linkingvardata.linkconss[2]);
}

TEST_F(GcgVarTest, RayMastervarIsRay) {
   MASTERVAR(var, vardata);
   vardata.data.mastervardata.isray = TRUE;
   ASSERT_EQ(TRUE, GCGmasterVarIsRay(&var));
}

TEST_F(GcgVarTest, NonRayMastervarIsNotRay) {
   MASTERVAR(var, vardata);
   vardata.data.mastervardata.isray = FALSE;
   ASSERT_EQ(FALSE, GCGmasterVarIsRay(&var));
}

TEST_F(GcgVarTest, MastervarGetOrigvars) {
   MASTERVAR(var, vardata);
   ORIGVAR(var2, vardata2);
   SCIP_VAR* vars[1] = { &var2 };
   vardata2.blocknr = 1;
   vardata.data.mastervardata.origvars = vars;
   vardata.data.mastervardata.norigvars = 1;
   ASSERT_EQ((SCIP_VAR**) &vars, GCGmasterVarGetOrigvars(&var));
}

TEST_F(GcgVarTest, MastervarGetNOrigvars) {
   MASTERVAR(var, vardata);
   ORIGVAR(var2, vardata2);
   SCIP_VAR* vars[1] = { &var2 };
   vardata.blocknr = 1;
   vardata.data.mastervardata.origvars = vars;
   vardata.data.mastervardata.norigvars = 0xDEAD;
   ASSERT_EQ(0xDEAD, GCGmasterVarGetNOrigvars(&var));
}

TEST_F(GcgVarTest, MastervarGetOrigvals) {
   MASTERVAR(var, vardata);
   vardata.data.mastervardata.origvals = (SCIP_Real*) 0xDEADBEEF;
   ASSERT_EQ((SCIP_Real*) 0xDEADBEEF, GCGmasterVarGetOrigvals(&var));
}

TEST_F(GcgVarTest, PricingvarGetOrigvars) {
   PRICINGVAR(var, vardata);
   vardata.data.pricingvardata.origvars = (SCIP_VAR**) 0xDEADBEEF;
   ASSERT_EQ((SCIP_VAR**) 0xDEADBEEF, GCGpricingVarGetOrigvars(&var));
}

TEST_F(GcgVarTest, PricingvarGetNOrigvars) {
   PRICINGVAR(var, vardata);
   vardata.data.pricingvardata.norigvars = 0xDEAD;
   ASSERT_EQ(0xDEAD, GCGpricingVarGetNOrigvars(&var));
}

TEST_F(GcgVarTest, LinkingVarInBlockIsInBlock) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkingvardata;
   SCIP_VAR* pricingvars[3] = {&lvar, NULL, &lvar};
   linkingvardata.pricingvars = pricingvars;
   linkingvardata.nblocks = 2;
   lvardata.data.origvardata.linkingvardata = &linkingvardata;
   ASSERT_EQ(TRUE, GCGisLinkingVarInBlock(&lvar, 2));
}

TEST_F(GcgVarTest, LinkingVarNotInBlockIsNotInBlock) {
   LINKINGVAR(lvar, lvardata);
   GCG_LINKINGVARDATA linkingvardata;
   SCIP_VAR* pricingvars[3] = {&lvar, NULL, &lvar};
   linkingvardata.pricingvars = pricingvars;
   linkingvardata.nblocks = 2;
   lvardata.data.origvardata.linkingvardata = &linkingvardata;
   ASSERT_EQ(FALSE, GCGisLinkingVarInBlock(&lvar, 1));
}

TEST_F(GcgVarTest, OriginalVarAddMasterVarWithReallocation)
{
   ORIGVAR(ovar, ovardata);
   MASTERVAR(mvar, mvardata);
   ovardata.data.origvardata.maxmastervars = 1;
   ovardata.data.origvardata.nmastervars = 1;
   SCIP_CALL_EXPECT(SCIPincludeRelaxGcg(scip));

   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.mastervars, 1));
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.mastervals, 1));
   ovardata.data.origvardata.mastervars[0] = (SCIP_VAR*) 0xDEADBEEF;
   ovardata.data.origvardata.mastervals[0] = 1.0;
   SCIP_CALL_EXPECT(GCGoriginalVarAddMasterVar(scip, &ovar, &mvar, 2.0));

   ASSERT_EQ((SCIP_VAR*) 0xDEADBEEF, ovardata.data.origvardata.mastervars[0]);
   ASSERT_EQ(1.0, ovardata.data.origvardata.mastervals[0]);
   ASSERT_EQ(&mvar, ovardata.data.origvardata.mastervars[1]);
   ASSERT_EQ(2.0, ovardata.data.origvardata.mastervals[1]);
   ASSERT_EQ(2, ovardata.data.origvardata.nmastervars);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.mastervars, ovardata.data.origvardata.maxmastervars);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.mastervals, ovardata.data.origvardata.maxmastervars);
}

TEST_F(GcgVarTest, OriginalVarAddMasterVarWithoutReallocation)
{
   ORIGVAR(ovar, ovardata);
   MASTERVAR(mvar, mvardata);
   ovardata.data.origvardata.maxmastervars = 2;
   ovardata.data.origvardata.nmastervars = 1;
   SCIP_CALL_EXPECT(SCIPincludeRelaxGcg(scip));
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.mastervars, 2));
   SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata.data.origvardata.mastervals, 2));
   ovardata.data.origvardata.mastervars[0] = (SCIP_VAR*) 0xDEADBEEF;
   ovardata.data.origvardata.mastervals[0] = 1.0;
   ovardata.data.origvardata.mastervars[1] = (SCIP_VAR*) 0xDEADBEEF;
   ovardata.data.origvardata.mastervals[1] = 2.0;
   SCIP_CALL_EXPECT(GCGoriginalVarAddMasterVar(scip, &ovar, &mvar, 2.0));

   ASSERT_EQ((SCIP_VAR*) 0xDEADBEEF, ovardata.data.origvardata.mastervars[0]);
   ASSERT_EQ(1.0, ovardata.data.origvardata.mastervals[0]);
   ASSERT_EQ(&mvar, ovardata.data.origvardata.mastervars[1]);
   ASSERT_EQ(2.0, ovardata.data.origvardata.mastervals[1]);
   ASSERT_EQ(2, ovardata.data.origvardata.nmastervars);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.mastervars, 2);
   SCIPfreeBlockMemoryArray(scip, &ovardata.data.origvardata.mastervals, 2);
}

TEST_F(GcgVarTest, OriginalVarRemoveExistingMasterVar)
{
   ORIGVAR(ovar, ovardata);
   MASTERVAR(mvar, mvardata);

   SCIP_VAR* vars[4] = {(SCIP_VAR*) 0xDEADBEEF, &mvar, (SCIP_VAR*) 0xDEADCAFF, NULL};
   SCIP_Real vals[4] = {1.0, 2.0, 3.0, -1.0};

   ovardata.data.origvardata.maxmastervars = 4;
   ovardata.data.origvardata.nmastervars = 3;
   ovardata.data.origvardata.mastervars = vars;
   ovardata.data.origvardata.mastervals = vals;

   SCIP_CALL_EXPECT(GCGoriginalVarRemoveMasterVar(scip, &ovar, &mvar));

   ASSERT_EQ((SCIP_VAR*) 0xDEADBEEF, ovardata.data.origvardata.mastervars[0]);
   ASSERT_EQ(1.0, ovardata.data.origvardata.mastervals[0]);
   ASSERT_EQ((SCIP_VAR*) 0xDEADCAFF, ovardata.data.origvardata.mastervars[1]);
   ASSERT_EQ(3.0, ovardata.data.origvardata.mastervals[1]);
   ASSERT_EQ(2, ovardata.data.origvardata.nmastervars);

}

TEST_F(GcgVarTest, OriginalVarRemoveNonExistingMasterVar)
{
   ORIGVAR(ovar, ovardata);
   MASTERVAR(mvar, mvardata);
   SCIP_VAR* vars[4] = {(SCIP_VAR*) 0xDEADBEEF, (SCIP_VAR*) 0xDEADCAFF, NULL, NULL};
   SCIP_Real vals[4] = {1.0, 3.0, -1.0, -1.0};
   ovardata.data.origvardata.maxmastervars = 4;
   ovardata.data.origvardata.nmastervars = 2;
   ovardata.data.origvardata.mastervars = vars;
   ovardata.data.origvardata.mastervals = vals;

   SCIP_CALL_EXPECT(GCGoriginalVarRemoveMasterVar(scip, &ovar, &mvar));

   ASSERT_EQ((SCIP_VAR*) 0xDEADBEEF, ovardata.data.origvardata.mastervars[0]);
   ASSERT_EQ(1.0, ovardata.data.origvardata.mastervals[0]);
   ASSERT_EQ((SCIP_VAR*) 0xDEADCAFF, ovardata.data.origvardata.mastervars[1]);
   ASSERT_EQ(3.0, ovardata.data.origvardata.mastervals[1]);
   ASSERT_EQ(2, ovardata.data.origvardata.nmastervars);
}

TEST_F(GcgVarTest, OriginalVarCreatePricingVar)
{
   SCIP_VAR* pricingvar = NULL;
   SCIP_VAR* ovar;
   SCIP_VARDATA ovardata;
   ovardata.blocknr = 0;
   ovardata.vartype = GCG_VARTYPE_ORIGINAL;
   ovardata.data.origvardata.linkingvardata = NULL;
   ovardata.data.origvardata.pricingvar = NULL;
   SCIP_CALL_EXPECT(SCIPcreateProbBasic(scip, "temp"));
   SCIP_CALL_EXPECT(SCIPcreateVarBasic(scip, &ovar, "test", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY));
   SCIPvarSetData(ovar, &ovardata);

   SCIP_CALL_EXPECT(GCGoriginalVarCreatePricingVar(scip, ovar, &pricingvar));

   ASSERT_NE((SCIP_VAR*) NULL, pricingvar);
   ASSERT_EQ(1, GCGvarIsPricing(pricingvar));
   ASSERT_EQ(0.0, SCIPvarGetLbGlobal(pricingvar));
   ASSERT_EQ(1.0, SCIPvarGetUbGlobal(pricingvar));
   ASSERT_EQ(0.0, SCIPvarGetObj(pricingvar));
   ASSERT_EQ(SCIP_VARTYPE_BINARY, SCIPvarGetType(pricingvar));
   ASSERT_EQ(1, GCGpricingVarGetNOrigvars(pricingvar));
   ASSERT_EQ(ovar, GCGpricingVarGetOriginalVar(pricingvar));

   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &pricingvar));
   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &ovar));
}

TEST_F(GcgVarTest, LinkingVarCreatePricingVar)
{
   SCIP_VAR* pricingvar = NULL;
   SCIP_VAR* ovar;
   SCIP_VARDATA ovardata;
   SCIP_CONS* linkcons = NULL;
   ovardata.blocknr = -2;
   ovardata.vartype = GCG_VARTYPE_ORIGINAL;
   ovardata.data.origvardata.pricingvar = NULL;

   SCIP_CALL_EXPECT(SCIPincludeConshdlrLinear(scip));
   SCIP_CALL_EXPECT(SCIPcreateProbBasic(scip, "temp"));
   SCIP_CALL_EXPECT(SCIPcreateVarBasic(scip, &ovar, "test", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY));
   SCIPvarSetData(ovar, &ovardata);

   SCIP_CALL_EXPECT(GCGlinkingVarCreatePricingVar(scip, 0, ovar, &pricingvar));
   GCGlinkingVarCreateMasterCons(scip,0, ovar, &linkcons);
   ASSERT_NE((SCIP_VAR*) NULL, pricingvar);
   ASSERT_NE((SCIP_CONS*) NULL, linkcons);
   ASSERT_EQ(0, SCIPgetNVarsLinear(scip, linkcons));

   ASSERT_EQ(1, GCGvarIsPricing(pricingvar));
   ASSERT_EQ(0.0, SCIPvarGetLbGlobal(pricingvar));
   ASSERT_EQ(1.0, SCIPvarGetUbGlobal(pricingvar));
   ASSERT_EQ(0.0, SCIPvarGetObj(pricingvar));
   ASSERT_EQ(SCIP_VARTYPE_BINARY, SCIPvarGetType(pricingvar));
   ASSERT_EQ(1, GCGpricingVarGetNOrigvars(pricingvar));
   ASSERT_EQ(ovar, GCGpricingVarGetOrigvars(pricingvar)[0]);

   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &pricingvar));
   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &ovar));
   SCIP_CALL_EXPECT(SCIPreleaseCons(scip, &linkcons));
}

TEST_F(GcgVarTest, CreateMasterVar)
{
   ORIGVAR(ovar1, ovard1);
   ORIGVAR(ovar2, ovard2);
   SCIP_VAR *newvar;
   SCIP_VAR *solvars[2];
   SCIP_Real solvals[2] = {2.0, -3.0};

   SCIP_VARDATA pvardata[2];
   SCIP_VARDATA* ovardata[2] = {&ovard1, &ovard2};
   SCIP_VAR* ovars[2] ={&ovar1, &ovar2};
   SCIP_VARDATA* mvardata;

   for( int i = 0; i < 2; ++i)
   {
      pvardata[i].blocknr = 0;
      pvardata[i].vartype = GCG_VARTYPE_PRICING;

      SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &pvardata[i].data.pricingvardata.origvars, 1));

      pvardata[i].data.pricingvardata.norigvars = 1;
      pvardata[i].data.pricingvardata.origvars[0] = ovars[i];
      ovardata[i]->data.origvardata.maxmastervars = 1;
      ovardata[i]->data.origvardata.nmastervars = 0;
      SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata[i]->data.origvardata.mastervars, 1));
      SCIP_CALL_EXPECT(SCIPallocBlockMemoryArray(scip, &ovardata[i]->data.origvardata.mastervals, 1));
   }


   SCIP_CALL_EXPECT(SCIPcreateProbBasic(scip, "temp"));
   SCIP_CALL_EXPECT(SCIPcreateVarBasic(scip, &(solvars[0]), "test", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY));
   SCIP_CALL_EXPECT(SCIPcreateVarBasic(scip, &(solvars[1]), "test2", -2.0, -1.0, -3.0, SCIP_VARTYPE_CONTINUOUS));
   SCIPvarSetData(solvars[0], &(pvardata[0]));
   SCIPvarSetData(solvars[1], &(pvardata[1]));
   SCIP_CALL_EXPECT(SCIPincludeRelaxGcg(scip));

   SCIP_CALL_EXPECT(GCGcreateMasterVar(scip, scip, scip, &newvar, "newname", 1.0, SCIP_VARTYPE_INTEGER, FALSE, 0, 2, solvals, solvars, DEC_DECMODE_DANTZIGWOLFE));

   ASSERT_NE((SCIP_VAR*)NULL, newvar);
   ASSERT_EQ(1.0, SCIPvarGetObj(newvar));
   ASSERT_EQ(SCIP_VARTYPE_INTEGER, SCIPvarGetType(newvar));
   ASSERT_EQ(TRUE, GCGvarIsMaster(newvar));
   ASSERT_EQ(2, GCGmasterVarGetNOrigvars(newvar));
   ASSERT_EQ(0, GCGvarGetBlock(newvar));
   ASSERT_EQ(FALSE, GCGmasterVarIsRay(newvar));
   for(int i = 0; i < 2; ++i)
   {
      ASSERT_EQ(ovars[i], GCGmasterVarGetOrigvars(newvar)[i]);
      ASSERT_EQ(solvals[i], GCGmasterVarGetOrigvals(newvar)[i]);
      ASSERT_EQ(1, GCGoriginalVarGetNMastervars(ovars[i]));
      ASSERT_EQ(newvar, GCGoriginalVarGetMastervars(ovars[i])[0]);
   }

   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &(solvars[0])));
   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &(solvars[1])));
   mvardata = SCIPvarGetData(newvar);
   SCIPfreeBlockMemoryArrayNull(scip, &mvardata->data.mastervardata.origvals, mvardata->data.mastervardata.norigvars);
   SCIPfreeBlockMemoryArrayNull(scip, &mvardata->data.mastervardata.origvars, mvardata->data.mastervardata.norigvars);
   SCIPfreeBlockMemory(scip, &mvardata);

   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &newvar));

   for(int i = 0; i < 2; ++i)
   {
      SCIPfreeBlockMemoryArray(scip, &pvardata[i].data.pricingvardata.origvars, pvardata[i].data.pricingvardata.norigvars);
      SCIPfreeBlockMemoryArray(scip, &ovardata[i]->data.origvardata.mastervals, 1);
      SCIPfreeBlockMemoryArray(scip, &ovardata[i]->data.origvardata.mastervars, 1);
   }
}

TEST_F(GcgVarTest, CreateInitialLinkingMasterVar)
{
   SCIP_VAR* mvar = NULL;
   SCIP_VAR* ovar;
   SCIP_VARDATA ovardata;
   SCIP_VARDATA *mvardata;
   ovardata.blocknr = -2;
   ovardata.vartype = GCG_VARTYPE_ORIGINAL;
   ovardata.data.origvardata.pricingvar = NULL;

   SCIP_CALL_EXPECT(SCIPcreateProbBasic(scip, "temp"));
   SCIP_CALL_EXPECT(SCIPcreateVarBasic(scip, &ovar, "test", 0.0, 1.0, 2.0, SCIP_VARTYPE_BINARY));
   SCIPvarSetData(ovar, &ovardata);
   SCIP_CALL_EXPECT(GCGcreateInitialMasterVar(scip, ovar, &mvar) );
   ASSERT_NE((SCIP_VAR*) NULL, mvar);
   ASSERT_EQ(0.0, SCIPvarGetLbGlobal(mvar));
   ASSERT_EQ(1.0, SCIPvarGetUbGlobal(mvar));
   ASSERT_EQ(2.0, SCIPvarGetObj(mvar));
   ASSERT_EQ(SCIP_VARTYPE_BINARY, SCIPvarGetType(mvar));

   ASSERT_EQ(TRUE, GCGvarIsMaster(mvar));
   ASSERT_EQ(FALSE, GCGmasterVarIsRay(mvar));
   ASSERT_EQ(ovar, GCGmasterVarGetOrigvars(mvar)[0]);
   ASSERT_EQ(1, GCGmasterVarGetNOrigvars(mvar));
   ASSERT_EQ(1.0, GCGmasterVarGetOrigvals(mvar)[0]);

   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &ovar));
   mvardata = SCIPvarGetData(mvar);
   SCIPfreeBlockMemoryArrayNull(scip, &mvardata->data.mastervardata.origvals, mvardata->data.mastervardata.norigvars);
   SCIPfreeBlockMemoryArrayNull(scip, &mvardata->data.mastervardata.origvars, mvardata->data.mastervardata.norigvars);
   SCIPfreeBlockMemory(scip, &mvardata);
   SCIP_CALL_EXPECT(SCIPreleaseVar(scip, &mvar));
}

TEST_F(GcgVarTest, SetCreationNode)
{
   MASTERVAR(var, vardata);
   vardata.creationnode = -1L;
   GCGsetCreationNode(scip, &var, 1L);
   ASSERT_EQ(1, vardata.creationnode);
}

TEST_F(GcgVarTest, GetCreationNode)
{
   MASTERVAR(var, vardata);
   vardata.creationnode = 1L;
   ASSERT_EQ(1L, GCGgetCreationNode(scip, &var));
}

TEST_F(GcgVarTest, SetCreationTime)
{
   MASTERVAR(var, vardata);
   vardata.creationtime = 0.0;
   GCGsetCreationTime(scip, &var, 1.0);
   ASSERT_EQ(1.0, vardata.creationtime);
}

TEST_F(GcgVarTest, GetCreationTime)
{
   MASTERVAR(var, vardata);
   vardata.creationtime = 1.0;
   ASSERT_EQ(1.0, GCGgetCreationTime(scip, &var));
}

TEST_F(GcgVarTest, SetIteration)
{
   MASTERVAR(var, vardata);
   vardata.iteration = -1;
   GCGsetIteration(scip, &var, 1L);
   ASSERT_EQ(1, vardata.iteration);
}

TEST_F(GcgVarTest, GetIteration)
{
   MASTERVAR(var, vardata);
   vardata.iteration = 1L;
   ASSERT_EQ(1, GCGgetIteration(scip, &var));
}

TEST_F(GcgVarTest, SetGap)
{
   MASTERVAR(var, vardata);
   vardata.gap = 0.0;
   GCGsetGap(scip, &var, 1.0);
   ASSERT_EQ(1.0, vardata.gap);
}

TEST_F(GcgVarTest, GetGap)
{
   MASTERVAR(var, vardata);
   vardata.gap = 1.0;
   ASSERT_EQ(1.0, GCGgetGap(scip, &var));
}

TEST_F(GcgVarTest, SetRedcost)
{
   MASTERVAR(var, vardata);
   vardata.redcost = 0.0;
   GCGsetRedcost(scip, &var, -1.0);
   ASSERT_EQ(-1.0, vardata.redcost);
}

TEST_F(GcgVarTest, GetRedcost)
{
   MASTERVAR(var, vardata);
   vardata.redcost = -1.0;
   ASSERT_EQ(-1.0, GCGgetRedcost(scip, &var));
}
