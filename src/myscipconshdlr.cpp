/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_StrongTriangles.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for StrongTriangles constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "myscipconshdlr.hpp"

bool ConshdlrStrongTriangles::is_constr_consistent_integer(SCIP *scip, SCIP_SOL *sol, const int i, const int j, const int k)
{
   auto ii = min(min(i, j), k);
   auto kk = max(max(i, j), k);
   auto jj = ii == i ? (min(j, k)) : (ii == j ? min(i, k) : min(i, j));
   auto uu = (*subgraph)[ii];
   auto vv = (*subgraph)[jj];
   auto ww = (*subgraph)[kk];
   // u_int8_t cntdntcare = 0;
   // if (instance->CUV[uu][vv] == instance->CUV[vv][uu])cntdntcare++;
   // if (instance->CUV[uu][ww] == instance->CUV[ww][uu])cntdntcare++;
   // if (instance->CUV[vv][ww] == instance->CUV[ww][vv])cntdntcare++;
   // if (cntdntcare > 2) return true;
   double lb = 0;
   double ub = 1;
   double val = 0;
   if (instance->pomatrix[uu][vv] == 0)
   {
      assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][jj]])));
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][jj]]);
   }
   else if (instance->pomatrix[uu][vv] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[vv][ww] == 0)
   {
      assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[jj][kk]])));
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[jj][kk]]);
   }
   else if (instance->pomatrix[vv][ww] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[uu][ww] == 0)
   {
      assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][kk]])));
      val -= SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][kk]]);
   }
   else if (instance->pomatrix[uu][ww] == 1)
   {
      lb += 1;
      ub += 1;
   }
   return lb <= val && ub >= val;
}

bool ConshdlrStrongTriangles::separate_if_efficatious(SCIP *scip, SCIP_SOL *sol, const int i, const int j, const int k, unsigned int enforce, SCIP_CONSHDLR *hdlr, double minviolation, SCIP_Result *result)
{
   auto ii = min(min(i, j), k);
   auto kk = max(max(i, j), k);
   auto jj = ii == i ? (min(j, k)) : (ii == j ? min(i, k) : min(i, j));
   auto uu = (*subgraph)[ii];
   auto vv = (*subgraph)[jj];
   auto ww = (*subgraph)[kk];
   if (stage == 0)
   {
      u_int8_t cntdntcare = 0;
      if (instance->CUV[uu][vv] == instance->CUV[vv][uu])
         cntdntcare++;
      if (instance->CUV[uu][ww] == instance->CUV[ww][uu])
         cntdntcare++;
      if (instance->CUV[vv][ww] == instance->CUV[ww][vv])
         cntdntcare++;
      if (cntdntcare > 2)
         return false;
   }
   double lb = 0;
   double ub = 1;
   double val = 0;
   vector<SCIP_VAR *> indices;
   vector<double> values;
   if (instance->pomatrix[uu][vv] == 0)
   {
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][jj]]);
      values.push_back(1.0);
      indices.push_back((*scipvars)[(*var_index)[ii][jj]]);
   }
   else if (instance->pomatrix[uu][vv] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[vv][ww] == 0)
   {
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[jj][kk]]);
      values.push_back(1.0);
      indices.push_back((*scipvars)[(*var_index)[jj][kk]]);
   }
   else if (instance->pomatrix[vv][ww] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[uu][ww] == 0)
   {
      val -= SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][kk]]);
      values.push_back(-1.0);
      indices.push_back((*scipvars)[(*var_index)[ii][kk]]);
   }
   else if (instance->pomatrix[uu][ww] == 1)
   {
      lb += 1;
      ub += 1;
   }
   if (lb - minviolation > val || ub + minviolation < val)
   {
      assert(values.size() > 0);
      SCIP_ROW *row;
      SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, hdlr, "", lb, ub, false, false, true));

      for (u_int32_t i = 0; i < values.size(); i++)
      {
         SCIPaddVarToRow(scip, row, indices[i], values[i]);
      }
      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      if (enforce || (SCIPisCutEfficacious(scip, sol, row) && SCIPisCutNew(scip, row)));
      {
         SCIP_Bool infeasible;
         SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));
         if (infeasible)
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
      }
      SCIP_CALL(SCIPreleaseRow(scip, &row));
      return true;
   }
   return false;
}

bool ConshdlrStrongTriangles::storecut(SCIP *scip, SCIP_SOL *sol, const int i, const int j, const int k, unsigned int enforce, SCIP_CONSHDLR *hdlr, double minviolation, map<double, vector<SCIP_ROW *>> &cuts)
{
   auto ii = min(min(i, j), k);
   auto kk = max(max(i, j), k);
   auto jj = ii == i ? (min(j, k)) : (ii == j ? min(i, k) : min(i, j));
   auto uu = (*subgraph)[ii];
   auto vv = (*subgraph)[jj];
   auto ww = (*subgraph)[kk];

   double lb = 0;
   double ub = 1;
   double val = 0;
   vector<SCIP_VAR *> indices;
   vector<double> values;
   if (instance->pomatrix[uu][vv] == 0)
   {
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][jj]]);
      values.push_back(1.0);
      indices.push_back((*scipvars)[(*var_index)[ii][jj]]);
   }
   else if (instance->pomatrix[uu][vv] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[vv][ww] == 0)
   {
      val += SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[jj][kk]]);
      values.push_back(1.0);
      indices.push_back((*scipvars)[(*var_index)[jj][kk]]);
   }
   else if (instance->pomatrix[vv][ww] == 1)
   {
      lb -= 1;
      ub -= 1;
   }
   if (instance->pomatrix[uu][ww] == 0)
   {
      val -= SCIPgetSolVal(scip, sol, (*scipvars)[(*var_index)[ii][kk]]);
      values.push_back(-1.0);
      indices.push_back((*scipvars)[(*var_index)[ii][kk]]);
   }
   else if (instance->pomatrix[uu][ww] == 1)
   {
      lb += 1;
      ub += 1;
   }
   if (lb - minviolation > val || ub + minviolation < val)
   {
      double violation = max(lb - val, val - ub);
      assert(values.size() > 0);
      SCIP_ROW *row;
      SCIP_CALL(SCIPcreateEmptyRowConshdlr(scip, &row, hdlr, "", lb, ub, false, false, true));

      for (u_int32_t i = 0; i < values.size(); i++)
      {
         SCIPaddVarToRow(scip, row, indices[i], values[i]);
      }
      SCIP_CALL(SCIPflushRowExtensions(scip, row));
      if (enforce || SCIPisCutEfficacious(scip, sol, row))
      {
         cuts[-violation].push_back(row);
      }
      else {
         SCIP_CALL(SCIPreleaseRow(scip, &row));
      }      
      return true;
   }
   return false;
}

bool ConshdlrStrongTriangles::is_consistent_integer(SCIP *scip, SCIP_SOL *sol)
{
   // vector<tuple<u_int32_t, int, u_int32_t>> eventqueue;
   // for (auto u : *subgraph)
   // {
   //    if (instance->al[u].size() == 1)
   //    {
   //       eventqueue.push_back({instance->al[u][instance->al[u].size() - 1] * 2 + 1, 0, (*index)[u]});
   //       eventqueue.push_back({instance->al[u][0] * 2, 1, (*index)[u]});
   //    }
   //    else
   //    {
   //       eventqueue.push_back({instance->al[u][instance->al[u].size() - 1] * 2, 0, (*index)[u]});
   //       eventqueue.push_back({instance->al[u][0] * 2, 1, (*index)[u]});
   //    }
   // }
   // sort(eventqueue.begin(), eventqueue.end());
   // set<u_int32_t> active;
   // for (auto event : eventqueue)
   // {
   //    auto k = get<2>(event);

   //    if (get<1>(event) == 0)
   //    {
   //       active.erase(k);
   //    }
   //    else
   //    {
   //       vector<u_int32_t> activevec(active.begin(), active.end());
   //       for (u_int32_t i = 0; i < activevec.size(); i++)
   //       {
   //          for (u_int32_t j = i + 1; j < activevec.size(); j++)
   //          {
   //             auto u = (*subgraph)[i];
   //             auto v = (*subgraph)[j];
   //             auto w = (*subgraph)[k];
   //             // if ((instance->CUV[u][v] != instance->CUV[v][u] || instance->CUV[v][w] != instance->CUV[w][v] || instance->CUV[u][w] != instance->CUV[w][u]))
   //             // {
   //                if (!this->is_constr_consistent_integer(scip, sol, i, j, k))
   //                {
   //                   return false;
   //                }
   //             // }
   //          }
   //       }
   //       active.insert(k);
   //    }
   // }

   for (u_int32_t i = 0; i < subgraph->size(); i++)
   {
      auto u = (*subgraph)[i];
      for (u_int32_t x = 0; x < intersecting[i].size(); x++)
      {

         for (u_int32_t y = x + 1; y < intersecting[i].size(); y++)
         {
            auto j = intersecting[i][x];
            auto k = intersecting[i][y];
            auto w = (*subgraph)[k];
            auto v = (*subgraph)[j];
            if (isintersect[j][k])
            {
               if (i < j && j < k)
               {
                  if (!this->is_constr_consistent_integer(scip, sol, i, j, k))
                  {
                     return false;
                  }
               }
            }
            else
            {
               if (!this->is_constr_consistent_integer(scip, sol, i, j, k))
               {
                  return false;
               }
            }
         }
      }
   }
   return true;
}

SCIP_RETCODE ConshdlrStrongTriangles::sepaStrongTriangles(SCIP *scip, SCIP_SOL *sol, SCIP_Bool enforce, SCIP_RESULT *result, SCIP_Conshdlr *hdlr)
{
   *result = SCIP_DIDNOTFIND;

   double minviolation = SCIPfeastol(scip);
   int cutcntweak = 0;
   int cutpoolweak = 0;
   int cutpoolstrong = 0;
   vector<tuple<u_int32_t, int, u_int32_t>> eventqueue;
   for (auto u : *subgraph)
   {
      if (instance->al[u].size() == 1)
      {
         eventqueue.push_back({instance->al[u][instance->al[u].size() - 1] + 1, 0, (*index)[u]});
         eventqueue.push_back({instance->al[u][0], 1, (*index)[u]});
      }
      else
      {
         eventqueue.push_back({instance->al[u][instance->al[u].size() - 1], 0, (*index)[u]});
         eventqueue.push_back({instance->al[u][0], 1, (*index)[u]});
      }
   }
   sort(eventqueue.begin(), eventqueue.end());
   set<u_int32_t> active;
   int cutpoolweakweak = 0;
   map<double, vector<SCIP_Row *>> cuts;
TRYAGAINAFTERSTAGECHANGE:
   for (auto event : eventqueue)
   {
      auto k = get<2>(event);

      if (get<1>(event) == 0)
      {
         active.erase(k);
      }
      else
      {
         vector<u_int32_t> activevec(active.begin(), active.end());
         for (u_int32_t ii = 0; ii < activevec.size(); ii++)
         {
            for (u_int32_t jj = ii + 1; jj < activevec.size(); jj++)
            {
               auto i = activevec[ii];
               auto j = activevec[jj];
               cutpoolweakweak++;
               // this->storecut(scip, sol, i, j, k, enforce, hdlr, minviolation, cuts);
               if (this->separate_if_efficatious(scip, sol, i, j, k, enforce, hdlr, minviolation, result))
               {
                  cutcntweak++;
                  if (*result == SCIP_CUTOFF)
                     return SCIP_OKAY;
               }
            }
         }
         active.insert(k);
      }
   }

   // for (const auto &kv : cuts)
   // {
   //    for (auto row : kv.second)
   //    {
   //       if (cutcntweak < 1000000)
   //       {
   //          cutcntweak++;
   //          SCIP_Bool infeasible;
   //          SCIP_CALL(SCIPaddRow(scip, row, FALSE, &infeasible));
   //          if (infeasible)
   //          {
   //             *result = SCIP_CUTOFF;
   //             return SCIP_OKAY;
   //          }
   //          else
   //             *result = SCIP_SEPARATED;
   //       }
   //       SCIP_CALL(SCIPreleaseRow(scip, &row));
   //    }
   // }

   if (cutcntweak)
   {
      debug(stage);
      debug(cutpoolweakweak);
      debug(cutcntweak);
      return SCIP_OKAY;
   }
   else if (stage == 0)
   {
      stage = 1;
      goto TRYAGAINAFTERSTAGECHANGE;
   }

   int cutcntstrong = 0;
   for (u_int32_t i = 0; i < subgraph->size(); i++)
   {
      for (u_int32_t x = 0; x < intersecting[i].size(); x++)
      {
         for (u_int32_t y = x + 1; y < intersecting[i].size(); y++)
         {
            auto j = intersecting[i][x];
            auto k = intersecting[i][y];
            auto u = (*subgraph)[i];
            auto v = (*subgraph)[j];
            auto w = (*subgraph)[k];
            if (isintersect[j][k])
            {
               if (i < j && j < k)
               {
                  cutpoolweak++;
                  if (this->separate_if_efficatious(scip, sol, i, j, k, enforce, hdlr, minviolation, result))
                  {
                     cutcntweak++;
                     if (*result == SCIP_CUTOFF)
                        return SCIP_OKAY;
                  }
               }
            }
            else
            {
               cutpoolstrong++;
               if (this->separate_if_efficatious(scip, sol, i, j, k, enforce, hdlr, minviolation, result))
               {
                  cutcntstrong++;
                  if (*result == SCIP_CUTOFF)
                     return SCIP_OKAY;
               }
            }
         }
      }
   }
   debug(cutpoolweak, cutpoolstrong);
   debug(cutcntweak, cutcntstrong);
   return SCIP_OKAY;
}

// void ConshdlrStrongTriangles::propagatetriangle(SCIP *scip, u_int32_t i, u_int32_t j, u_int32_t k, SCIP_RESULT *result, SCIP_Bool *tightened, SCIP_Bool *infeasible, SCIP_CONSHDLR *hdlr)
// {
//    auto ii = min(min(i, j), k);
//    auto kk = max(max(i, j), k);
//    auto jj = ii == i ? (min(j, k)) : (ii == j ? min(i, k) : min(i, j));

//    auto ibj = SCIPvarGetLbLocal((*scipvars)[(*this->var_index)[ii][jj]]) > 0.5;
//    auto jbi = SCIPvarGetUbLocal((*scipvars)[(*this->var_index)[ii][jj]]) < 0.5;
//    auto ibk = SCIPvarGetLbLocal((*scipvars)[(*this->var_index)[ii][kk]]) > 0.5;
//    auto kbi = SCIPvarGetUbLocal((*scipvars)[(*this->var_index)[ii][kk]]) < 0.5;
//    auto jbk = SCIPvarGetLbLocal((*scipvars)[(*this->var_index)[jj][kk]]) > 0.5;
//    auto kbj = SCIPvarGetUbLocal((*scipvars)[(*this->var_index)[jj][kk]]) < 0.5;
//    if (ibj && jbk && !ibk)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[ii][kk]]), 1, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//       {
//          return;
//       }
//    }
//    if (ibk && kbj && !ibj)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[ii][jj]]), 1, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//          return;
//    }
//    if (jbi && ibk && !jbk)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[jj][kk]]), 1, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//          return;
//    }
//    if (kbj && jbi && !kbi)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[ii][kk]]), 0, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//          return;
//    }
//    if (jbk && kbi && !jbi)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[ii][jj]]), 0, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//          return;
//    }
//    if (kbi && ibj && !kbj)
//    {
//       SCIP_Cons *cons = 0;
//       SCIPinferBinvarCons(scip, ((*scipvars)[(*this->var_index)[jj][kk]]), 0, dummycons, 0, infeasible, tightened);
//       if (*infeasible)
//          return;
//    }
// }

// SCIP_DECL_CONSPROP(ConshdlrStrongTriangles::scip_prop)
// {
//    // debug("in propagation", nusefulconss, nmarkedconss, proptiming);
//    int i = 0;
//    for (auto var: *scipvars) {
//       i++;
//       if (SCIPvarGetLbLocal(var) > 0.5){
//          debug(i);
//          break;
//       }
//    }

//    *result = SCIP_DIDNOTFIND;
//    int ntightened = 0;
//    for (u_int32_t i = 0; i < subgraph->size(); i++)
//    {
//       for (u_int32_t x = 0; x < intersecting[i].size(); x++)
//       {
//          for (u_int32_t y = x + 1; y < intersecting[i].size(); y++)
//          {
//             auto j = intersecting[i][x];
//             auto k = intersecting[i][y];
//             auto u = (*subgraph)[i];
//             auto v = (*subgraph)[j];
//             auto w = (*subgraph)[k];
//             if (isintersect[j][k])
//             {
//                if (i < j && j < k)
//                {
//                   SCIP_Bool tightenend = false;
//                   SCIP_Bool infeasible = false;
//                   propagatetriangle(scip, i, j, k, result, &tightenend, &infeasible, conshdlr);
//                   if (infeasible)
//                   {
//                      *result = SCIP_CUTOFF;
//                      debug("cutoff");
//                      return SCIP_OKAY;
//                   }
//                   if (tightenend)
//                   {
//                      ntightened++;
//                   }
//                }
//             }
//             else
//             {
//                SCIP_Bool tightenend = false;
//                SCIP_Bool infeasible = false;
//                propagatetriangle(scip, i, j, k, result, &tightenend, &infeasible, conshdlr);
//                if (infeasible)
//                {
//                   *result = SCIP_CUTOFF;
//                   debug("cutoff");
//                   return SCIP_OKAY;
//                }
//                if (tightenend)
//                {
//                   ntightened++;
//                }
//             }
//          }
//       }
//    }
//    debug(ntightened);
//    if (ntightened > 0)
//       *result = SCIP_REDUCEDDOM;
//    return SCIP_OKAY;
// }

SCIP_DECL_CONSSEPALP(ConshdlrStrongTriangles::scip_sepalp)
{ /*lint --e{715}*/
   SCIP_CALL(this->sepaStrongTriangles(scip, NULL, false, result, conshdlr));
   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
SCIP_DECL_CONSSEPASOL(ConshdlrStrongTriangles::scip_sepasol)
{ /*lint --e{715}*/
   SCIP_CALL(this->sepaStrongTriangles(scip, sol, false, result, conshdlr));
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(ConshdlrStrongTriangles::scip_enfolp)
{ /*lint --e{715}*/
   SCIP_CALL(this->sepaStrongTriangles(scip, NULL, true, result, conshdlr));

   if (*result == SCIP_DIDNOTFIND)
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

// /** constraint enforcing method of constraint handler for relaxation solutions */
// SCIP_DECL_CONSENFORELAX(ConshdlrStrongTriangles::scip_enforelax)
// { /*lint --e{715}*/
//    debug("doing constraint enforcing for relaxation");
//    SCIPerrorMessage("method of StrongTriangles constraint handler not implemented yet\n");
//    SCIPABORT(); /*lint --e{527} --e{715}*/
//    return SCIP_OKAY;
// }

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(ConshdlrStrongTriangles::scip_enfops)
{ /*lint --e{715}*/
   SCIP_CALL(this->sepaStrongTriangles(scip, NULL, true, result, conshdlr));
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
SCIP_DECL_CONSCHECK(ConshdlrStrongTriangles::scip_check)
{ /*lint --e{715}*/
   if (is_consistent_integer(scip, sol))
   {
      *result = SCIP_FEASIBLE;
   }
   else
   {
      *result = SCIP_INFEASIBLE;
   }
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(ConshdlrStrongTriangles::scip_lock)
{ /*lint --e{715}*/

   for (auto var : *scipvars)
   {
      SCIP_CALL(SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg));
   }

   return SCIP_OKAY;
}

/** creates and captures a StrongTriangles constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsStrongTriangles(
    SCIP *scip,              /**< SCIP data structure */
    SCIP_CONS **cons,        /**< pointer to hold the created constraint */
    const char *name,        /**< name of constraint */
    int nvars,               /**< number of variables in the constraint */
    SCIP_VAR **vars,         /**< array with variables of constraint entries */
    SCIP_Real *coefs,        /**< array with coefficients of constraint entries */
    SCIP_Real lhs,           /**< left hand side of constraint */
    SCIP_Real rhs,           /**< right hand side of constraint */
    SCIP_Bool initial,       /**< should the LP relaxation of constraint be in the initial LP?
                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
    SCIP_Bool separate,      /**< should the constraint be separated during LP processing?
                              *   Usually set to TRUE. */
    SCIP_Bool enforce,       /**< should the constraint be enforced during node processing?
                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
    SCIP_Bool check,         /**< should the constraint be checked for feasibility?
                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
    SCIP_Bool propagate,     /**< should the constraint be propagated during node processing?
                              *   Usually set to TRUE. */
    SCIP_Bool local,         /**< is constraint only valid locally?
                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
    SCIP_Bool modifiable,    /**< is constraint modifiable (subject to column generation)?
                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                              *   adds coefficients to this constraint. */
    SCIP_Bool dynamic,       /**< is constraint subject to aging?
                              *   Usually set to FALSE. Set to TRUE for own cuts which
                              *   are separated as constraints. */
    SCIP_Bool removable,     /**< should the relaxation be removed from the LP due to aging or cleanup?
                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
    SCIP_Bool stickingatnode /**< should the constraint always be kept at the node where it was added, even
                              *   if it may be moved to a more global node?
                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
)
{
   SCIP_CONSHDLR *conshdlr;
   SCIP_CONSDATA *consdata;

   SCIPerrorMessage("method of StrongTriangles constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the StrongTriangles constraint handler */
   // TODO: this might be slow, check how often called
   conshdlr = SCIPfindConshdlr(scip, "Strong triangles constraint handler");
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("StrongTriangles constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                            local, modifiable, dynamic, removable, stickingatnode));

   return SCIP_OKAY;
}

// SCIP_DECL_CONSINITLP(ConshdlrStrongTriangles::init)
// {
//    SCIPcreateConsStrongTriangles(scip, &dummycons, "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
// }