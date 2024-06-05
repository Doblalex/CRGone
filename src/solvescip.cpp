#include "solvescip.hpp"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "myscipconshdlr.hpp"
#include "myscipheuristic.hpp"

#define SCIP_CALL_ERROR(x)                 \
    do                                     \
    {                                      \
        SCIP_RETCODE _restat_;             \
        if ((_restat_ = (x)) != SCIP_OKAY) \
        {                                  \
            SCIPprintError(_restat_);      \
            return -1;                     \
        }                                  \
    } while (FALSE)

void addAllRows(OLCMInstance &instance, SCIP *scip, vector<u_int32_t> &subgraph, vector<SCIP_VAR *> &scipvars, vector<vector<u_int32_t>> &var_index)
{
    vector<vector<bool>> isintersect = vector<vector<bool>>(subgraph.size(), vector<bool>(subgraph.size(), false));
    vector<vector<u_int32_t>> intersecting = vector<vector<u_int32_t>>(subgraph.size());
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        auto u = subgraph[i];

        for (u_int32_t j = 0; j < subgraph.size(); j++)
        {
            auto v = subgraph[j];
            auto sa = instance.al[u][0];
            auto ea = instance.al[u][instance.al[u].size() - 1];
            auto sb = instance.al[v][0];
            auto eb = instance.al[v][instance.al[v].size() - 1];
            if (sa == eb || sb == ea)
                continue;
            // isintersect[i][j] = true;
            if ((sa <= sb && sb <= ea) || (sb <= sa && sa <= eb))
            {
                isintersect[i][j] = true;
            }
        }
    }

    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        auto u = subgraph[i];

        for (u_int32_t j = 0; j < subgraph.size(); j++)
        {
            if (i != j)
            {
                if (isintersect[i][j])
                    intersecting[i].push_back(j);
            }
        }
    }
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        for (u_int32_t x = 0; x < intersecting[i].size(); x++)
        {
            for (u_int32_t y = x + 1; y < intersecting[i].size(); y++)
            {
                auto j = intersecting[i][x];
                auto k = intersecting[i][y];
                auto u = subgraph[i];
                auto v = subgraph[j];
                auto w = subgraph[k];
                if (isintersect[j][k])
                {
                    if (i < j && j < k)
                    {
                        SCIP_VAR *vars[] = {scipvars[var_index[i][j]], scipvars[var_index[j][k]], scipvars[var_index[i][k]]};
                        double obj[] = {1, 1, -1};
                        SCIP_CONS *cons;
                        SCIPcreateConsBasicLinear(scip,
                                                  &cons,
                                                  "",
                                                  3,
                                                  vars,
                                                  obj,
                                                  0,
                                                  1);
                        SCIPaddCons(scip, cons);
                        SCIPreleaseCons(scip, &cons);
                    }
                }
                else
                {
                }
            }
        }
    }
}

SCIP_Retcode whatever(OLCMInstance &instance, vector<u_int32_t> &subgraph, vector<u_int32_t> &multiplicities, vector<u_int32_t> &ans)
{
    SCIP *scip = nullptr;
    SCIP_CALL(SCIPcreate(&scip)); // Creating the SCIP environment

    /* include default plugins */
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
#ifndef MYLOCAL
    SCIPsetMessagehdlrQuiet(scip, 1);
#endif

    // Creating the SCIP Problem.
    SCIP_CALL(SCIPcreateProbBasic(scip, "olcm"));
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
    // SCIP_CALL(SCIPsetIntParam(scip, "separating/maxcutsroot", 2000));
    SCIP_CALL(SCIPsetIntParam(scip, "separating/maxstallroundsroot", -1));
    SCIP_CALL(SCIPsetIntParam(scip, "lp/threads", 1));
    SCIP_CALL(SCIPsetIntParam(scip, "parallel/maxnthreads", 1));
    // SCIP_CALL(SCIPsetIntParam(scip, "lp/rowagelimit", 0));
    // SCIP_CALL(SCIPsetSubscipsOff(scip, TRUE));
    // SCIP_CALL( SCIPsetIntParam(scip, "separating/cutagelimit", 0) );
    // SCIP_CALL( SCIPsetIntParam(scip, "constraints/agelimit", 1) );
    debug("Solving ILP witn number of nodes", subgraph.size());

    static vector<u_int32_t> index(instance.nL2);
    vector<vector<u_int32_t>> var_index(subgraph.size(), vector<u_int32_t>(subgraph.size(), 0));
    vector<pair<u_int32_t, u_int32_t>> var_vals;
    vector<SCIP_VAR *> scipvars;

    for (u_int32_t i = 0; i < subgraph.size(); i++)
        index[subgraph[i]] = i;

    int32_t lb = 0;
    int cnt = 0;
    int cntuvvusame = 0;
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        for (u_int32_t j = i + 1; j < subgraph.size(); j++)
        {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (instance.pomatrix[u][v] == 0)
            {
                if (instance.CUV[u][v] == instance.CUV[v][u])
                    cntuvvusame++;
                // we do not know the order
                var_index[i][j] = cnt++;
                var_vals.push_back({i, j});
                int32_t mult = ((int32_t)multiplicities[i]) * ((int32_t)multiplicities[j]);
                int32_t obj = mult * (((int32_t)instance.CUV[u][v]) - ((int32_t)instance.CUV[v][u]));
                lb += ((int32_t)instance.CUV[v][u]) * mult;
                string s = to_string(i) + "_" + to_string(j);
                SCIP_VAR *x = nullptr;

                SCIP_CALL(SCIPcreateVarBasic(scip,
                                             &x,                    // returns new index
                                             "",                    // name
                                             0.0,                   // lower bound
                                             1.0,                   // upper bound
                                             obj,                   // objective
                                             SCIP_VARTYPE_BINARY)); // variable type
                SCIP_CALL(SCIPaddVar(scip, x));
                scipvars.push_back(x);
            }
        }
    }
    debug(cntuvvusame, scipvars.size());

    string collector;

#ifdef MYLOCAL
#else
    // TODO No logging
    BufferStdout swapstdout(collector);
#endif

    int n = subgraph.size();
    int nchoose3 = (n * (n - 1) * (n - 2)) / 6;
    debug(nchoose3);
    // addAllRows(instance, scip, subgraph, scipvars, var_index);
    SCIP_CALL(SCIPincludeObjConshdlr(scip, new ConshdlrStrongTriangles(scip, &var_index, &subgraph, &instance, &index, &scipvars), true));
    SCIP_CALL(SCIPincludeObjHeur(scip, new MyScipHeur(scip, &var_index, &subgraph, &instance, &index, &scipvars, &var_vals), true));
    cerr << "solving" << endl;
    // Solving the problem
    SCIP_CALL(SCIPsolve(scip));
    cerr << "solved" << endl;

    SCIP_SOL *sol = nullptr;

    sol = SCIPgetBestSol(scip);

    instance.LB += lb + SCIPgetSolOrigObj(scip, sol);

    // std::cout << "x1: " << SCIPgetSolVal(scip, sol, x1) << " " << "x2: " << SCIPgetSolVal(scip, sol, x2) << "\n";
    vector<vector<u_int32_t>> graph(subgraph.size());
    vector<u_int32_t> incoming(subgraph.size(), 0);
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        auto u = subgraph[i];
        for (u_int32_t j = i + 1; j < subgraph.size(); j++)
        {
            auto v = subgraph[j];
            if (instance.pomatrix[u][v] == 1)
            {
                graph[i].push_back(j);
                incoming[j]++;
            }
            else if (instance.pomatrix[u][v] == -1)
            {
                graph[j].push_back(i);
                incoming[i]++;
            }
            else if (instance.CUV[u][v] != instance.CUV[v][u])
            {
                auto val = SCIPgetSolVal(scip, sol, scipvars[var_index[i][j]]);
                if (val > 0.5)
                {
                    graph[i].push_back(j);
                    incoming[j]++;
                }
                else
                {
                    graph[j].push_back(i);
                    incoming[i]++;
                }
            }
        }
    }
    vector<u_int32_t> q;
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        if (incoming[i] == 0)
            q.push_back(i);
    }
    while (q.size())
    {
        auto i = q[q.size() - 1];
        q.pop_back();
        ans.push_back(subgraph[i]);
        for (auto j : graph[i])
        {
            if ((--incoming[j]) == 0)
            {
                q.push_back(j);
            }
        }
    }
    debug(ans.size(), subgraph.size());
    assert(ans.size() == subgraph.size());

    for (auto var : scipvars)
    {
        SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    SCIP_CALL(SCIPfree(&scip));

    return SCIP_OKAY;
}
vector<u_int32_t> solve_ilp(OLCMInstance &instance, vector<u_int32_t> &subgraph, vector<u_int32_t> &multiplicities)
{
    vector<u_int32_t> ans;
    if (subgraph.size() == 1)
    {
        ans.push_back(subgraph[0]);
        return ans;
    }
    whatever(instance, subgraph, multiplicities, ans);
    return ans;
}