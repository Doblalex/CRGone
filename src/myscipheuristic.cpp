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

/**@file   MyScipHeur.cpp
 * @brief  2-Optimum - combinatorial improvement heuristic for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <iostream>

#include "objscip/objscip.h"
#include "myscipheuristic.hpp"

using namespace std;

/** default constructor */
MyScipHeur::MyScipHeur(
    SCIP *scip,
    vector<vector<u_int32_t>> *var_index,
    vector<u_int32_t> *subgraph,
    OLCMInstance *instance,
    vector<u_int32_t> *index,
    vector<SCIP_Var *> *scipvars,
    vector<pair<u_int32_t, u_int32_t>> *var_vals)
    : ObjHeur(scip, "myrounding", "My rounding heuristic", 'K', 100, 1, 0, -1, SCIP_HEURTIMING_AFTERLPNODE | SCIP_HEURTIMING_DURINGLPLOOP, false)
{
    this->subgraph = subgraph;
    this->instance = instance;
    this->index = index;
    this->scipvars = scipvars;
    this->var_index = var_index;
    this->var_vals = var_vals;
    this->fixeddirection = vector<vector<bool>>(subgraph->size(), vector<bool>(subgraph->size(), false));
    vector<u_int32_t> incoming(subgraph->size(), 0);
    vector<vector<u_int32_t>> g(subgraph->size());
    sortoffset = vector<u_int32_t>(subgraph->size(), 0);
    for (u_int32_t i = 0; i < subgraph->size(); i++)
    {
        for (u_int32_t j = i + 1; j < subgraph->size(); j++)
        {
            auto val = instance->pomatrix[(*subgraph)[i]][(*subgraph)[j]];
            if (val == 1)
            {
                g[i].push_back(j);
                incoming[j]++;
                sortoffset[j]++;
            }
            else if (val == -1)
            {
                incoming[i]++;
                g[j].push_back(i);
                sortoffset[i]++;
            }
        }
    }

    vector<u_int32_t> sources;
    for (u_int32_t i = 0; i < subgraph->size(); i++)
    {
        if (incoming[i] == 0)
        {
            sources.push_back(i);
        }
    }
    while (sources.size())
    {
        vector<u_int32_t> newsources;
        for (auto u : sources)
        {
            someorder.push_back(u);
            for (auto v : g[u])
            {
                incoming[v]--;
                if (incoming[v] == 0)
                    newsources.push_back(v);
            }
        }
        sources = newsources;
    }
}

SCIP_DECL_HEURFREE(MyScipHeur::scip_free)
{ /*lint --e{715}*/
    return SCIP_OKAY;
}

SCIP_DECL_HEURINIT(MyScipHeur::scip_init)
{ /*lint --e{715}*/
    SCIP_CALL(SCIPcreateSol(scip, &sol, heur));
    return SCIP_OKAY;
}

SCIP_DECL_HEUREXIT(MyScipHeur::scip_exit)
{ /*lint --e{715}*/
    SCIP_CALL(SCIPfreeSol(scip, &sol));
    return SCIP_OKAY;
}

SCIP_DECL_HEURINITSOL(MyScipHeur::scip_initsol)
{
    return SCIP_OKAY;
}

SCIP_DECL_HEUREXITSOL(MyScipHeur::scip_exitsol)
{
    return SCIP_OKAY;
}

/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(MyScipHeur::scip_exec)
{
    SCIP_SOL *newsol;
    assert(result != NULL);
    /* since the timing is SCIP_HEURTIMING_AFTERLPNODE, the current node should have an LP */
    assert(SCIPhasCurrentNodeLP(scip));

    *result = SCIP_DIDNOTRUN;

    /* only call heuristic, if an optimal LP solution is at hand */
    if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL)
        return SCIP_OKAY;

    /* get the working solution from heuristic's local data */
    assert(sol != NULL);

    /* copy the current LP solution to the working solution */
    SCIP_CALL(SCIPlinkLPSol(scip, sol));

    *result = SCIP_DIDNOTFIND;

    /* allocate local memory */
    SCIP_CALL(SCIPcreateSol(scip, &newsol, heur));

    vector<double> sortby(subgraph->size(), 0.0);
    double solfixedvars = 0;
    for (u_int32_t i = 0; i < subgraph->size(); i++)
    {
        sortby[i] += sortoffset[i];
    }
    for (u_int32_t i = 0; i < var_vals->size(); i++)
    {
        auto p = (*var_vals)[i];
        auto var = (*scipvars)[i];
        sortby[p.first] += (1.0 - SCIPgetSolVal(scip, sol, var));
        sortby[p.second] += SCIPgetSolVal(scip, sol, var);
    }

    vector<u_int32_t> order_index(subgraph->size());
    vector<u_int32_t> inds(someorder.begin(), someorder.end());
    for (u_int32_t i = 1; i < subgraph->size(); i++)
    {
        auto j = i;
        while (j > 0 && sortby[inds[j]] <= sortby[inds[j - 1]] && instance->pomatrix[(*subgraph)[inds[j - 1]]][(*subgraph)[inds[j]]] != 1)
        {
            swap(inds[j - 1], inds[j]);
            j--;
        }
    }
    for (u_int32_t i = 0; i < subgraph->size(); i++)
    {
        order_index[inds[i]] = i;
    }

    double newsolutionvalue = 0;
    for (u_int32_t x = 0; x < var_vals->size(); x++)
    {
        auto p = (*var_vals)[x];
        auto i = p.first;
        auto j = p.second;
        auto u = (*subgraph)[i];
        assert(i != j);
        auto var = (*scipvars)[x];
        if (SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
        {
            if (SCIPvarGetLbLocal(var) > 0.5 && order_index[i] > order_index[j])
            {
                SCIP_CALL(SCIPfreeSol(scip, &newsol));
                return SCIP_OKAY;
            }
            if (SCIPvarGetUbLocal(var) < 0.5 && order_index[i] < order_index[j])
            {
                SCIP_CALL(SCIPfreeSol(scip, &newsol));
                return SCIP_OKAY;
            }
        }
        if (order_index[i] < order_index[j])
        {
            SCIPsetSolVal(scip, newsol, (*scipvars)[x], 1.0);
        }
        else
        {
            SCIPsetSolVal(scip, newsol, (*scipvars)[x], 0.0);
        }
    }
    SCIP_Bool success = false;

    /* due to construction we already know, that the solution will be feasible */
    SCIP_CALL(SCIPaddSol(scip, newsol, &success));
    // SCIP_CALL(SCIPtrySol(scip, newsol, FALSE, FALSE, FALSE, FALSE, FALSE, &success));
    if (success)
        *result = SCIP_FOUNDSOL;
    SCIP_CALL(SCIPfreeSol(scip, &newsol));
    return SCIP_OKAY;
}

// /** clone method which will be used to copy a objective plugin */
// SCIP_DECL_HEURCLONE(scip::ObjCloneable *MyScipHeur::clone) /*lint !e665*/
// {
//     return new MyScipHeur(scip);
// }
