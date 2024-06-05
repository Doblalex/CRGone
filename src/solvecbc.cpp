#include "solvecbc.hpp"
#include "mycbcheuristic.hpp"
#include "CoinPragma.hpp"
// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"
// For all different
#include "CbcBranchCut.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchAllDifferent.hpp"
#include "CbcCutGenerator.hpp"
#include "CglAllDifferent.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinTime.hpp"
#include "OsiAuxInfo.hpp"
#include "CglStored.hpp"
#include "CglPreProcess.hpp"
#include "CbcSolver.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicDive.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcStrategy.hpp"
#include "CbcFeasibilityBase.hpp"

typedef struct
{
    vector<vector<u_int32_t>> *var_index;
    vector<u_int32_t> *subgraph;
    OLCMInstance *instance;
    vector<u_int32_t> *index;
} CutInfo;

class CglTriangles : public CglCutGenerator
{
public:
    CglTriangles(const CutInfo *info);

    CglTriangles(const CglTriangles &rhs);

    virtual CglCutGenerator *clone() const;

    virtual void
    generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                 const CglTreeInfo info = CglTreeInfo());

    const CutInfo *info;
    vector<vector<bool>> isintersect;
    vector<vector<u_int32_t>> intersecting;
};

CglTriangles::CglTriangles(const CutInfo *i)
{
    this->info = i;
    auto subgraph = *this->info->subgraph;
    auto instance = this->info->instance;
    isintersect = vector<vector<bool>>(subgraph.size(), vector<bool>(subgraph.size(), false));
    intersecting = vector<vector<u_int32_t>>(subgraph.size());
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        auto u = subgraph[i];

        for (u_int32_t j = 0; j < subgraph.size(); j++)
        {
            auto v = subgraph[j];
            auto sa = instance->al[u][0];
            auto ea = instance->al[u][instance->al[u].size() - 1];
            auto sb = instance->al[v][0];
            auto eb = instance->al[v][instance->al[v].size() - 1];
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
}

CglTriangles::CglTriangles(const CglTriangles &rhs)
{
    this->info = rhs.info;
    this->isintersect = isintersect;
    this->intersecting = intersecting;
}

CglCutGenerator *CglTriangles::clone() const
{
    return new CglTriangles(this->info);
}

void CglTriangles::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
                                const CglTreeInfo tinfo)
{
    debug("adding cuts");
    const double *s = si.getColSolution();
    const auto var_index = *info->var_index;
    const auto subgraph = *info->subgraph;
    const auto instance = info->instance;
    const auto index = *info->index;
    int cutcntweak = 0;
    bool addedWeakCuts = false;

    static bool addedcuts = false;

    vector<tuple<u_int32_t, int, u_int32_t>> eventqueue;
    for (auto u : subgraph)
    {
        if (instance->al[u].size() == 1)
        {
            eventqueue.push_back({instance->al[u][instance->al[u].size() - 1] * 2 + 1, 0, index[u]});
            eventqueue.push_back({instance->al[u][0] * 2, 1, index[u]});
        }
        else
        {
            eventqueue.push_back({instance->al[u][instance->al[u].size() - 1] * 2, 0, index[u]});
            eventqueue.push_back({instance->al[u][0] * 2, 1, index[u]});
        }
    }
    sort(eventqueue.begin(), eventqueue.end());
    set<u_int32_t> active;
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
            for (u_int32_t i = 0; i < activevec.size(); i++)
            {
                for (u_int32_t j = i + 1; j < activevec.size(); j++)
                {
                    auto kk = k;
                    auto ii = activevec[i];
                    auto jj = activevec[j];
                    if (kk < ii)
                    {
                        swap(ii, kk);
                    }
                    if (kk < jj)
                    {
                        swap(jj, kk);
                    }
                    auto uu = subgraph[ii];
                    auto vv = subgraph[jj];
                    auto ww = subgraph[kk];
                    // if (instance->CUV[uu][vv] == instance->CUV[vv][uu] && instance->CUV[uu][ww] == instance->CUV[ww][uu] && instance->CUV[vv][ww] == instance->CUV[ww][vv])
                    if (instance->pomatrix[uu][vv] != 0)
                        continue;
                    if (instance->pomatrix[uu][ww] != 0)
                        continue;
                    if (instance->pomatrix[vv][ww] != 0)
                        continue;
                    vector<int> indices;
                    vector<double> values;
                    int lb = 0;
                    int ub = 1;
                    if (instance->pomatrix[uu][vv] == 0)
                    {
                        indices.push_back(var_index[ii][jj]);
                        values.push_back(1.0);
                    }
                    else if (instance->pomatrix[uu][vv] == 1)
                    {
                        lb -= 1;
                        ub -= 1;
                    }
                    if (instance->pomatrix[vv][ww] == 0)
                    {
                        indices.push_back(var_index[jj][kk]);
                        values.push_back(1.0);
                    }
                    else if (instance->pomatrix[vv][ww] == 1)
                    {
                        lb -= 1;
                        ub -= 1;
                    }
                    if (instance->pomatrix[uu][ww] == 0)
                    {
                        indices.push_back(var_index[ii][kk]);
                        values.push_back(-1.0);
                    }
                    else if (instance->pomatrix[uu][ww] == 1)
                    {
                        lb += 1;
                        ub += 1;
                    }
                    if (indices.size() == 0)
                        continue;
                    OsiRowCut cut;
                    cut.setLb(lb);
                    cut.setUb(ub);
                    cut.setRow(indices.size(), indices.data(), values.data());
                    if (cut.violated(s) > 1.0e-5)
                    {
                        cs.insert(cut);
                        addedWeakCuts = true;
                        cutcntweak++;
                    }
                }
            }
            active.insert(k);
        }
    }
    if (addedWeakCuts)
    {
        debug(cutcntweak);
        return;
    }

    int strongcutcnt = 0;
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        auto u = subgraph[i];
        for (u_int32_t x = 0; x < intersecting[i].size(); x++)
        {
            for (u_int32_t y = x + 1; y < intersecting[i].size(); y++)
            {
                auto j = intersecting[i][x];
                auto k = intersecting[i][y];
                auto v = subgraph[j];
                auto w = subgraph[k];

                if (isintersect[j][k])
                {
                    if (i < j && j < k) //&& (instance->CUV[u][v] != instance->CUV[v][u] || instance->CUV[v][w] != instance->CUV[w][v] || instance->CUV[u][w] != instance->CUV[w][u]))
                    {
                        vector<int> indices;
                        vector<double> values;
                        int lb = 0;
                        int ub = 1;
                        if (instance->pomatrix[u][v] == 0)
                        {
                            indices.push_back(var_index[i][j]);
                            values.push_back(1.0);
                        }
                        else if (instance->pomatrix[u][v] == 1)
                        {
                            lb -= 1;
                            ub -= 1;
                        }
                        if (instance->pomatrix[v][w] == 0)
                        {
                            indices.push_back(var_index[j][k]);
                            values.push_back(1.0);
                        }
                        else if (instance->pomatrix[v][w] == 1)
                        {
                            lb -= 1;
                            ub -= 1;
                        }
                        if (instance->pomatrix[u][w] == 0)
                        {
                            indices.push_back(var_index[i][k]);
                            values.push_back(-1.0);
                        }
                        else if (instance->pomatrix[u][w] == 1)
                        {
                            lb += 1;
                            ub += 1;
                        }
                        if (indices.size() == 0)
                            continue;
                        OsiRowCut cut;
                        cut.setLb(lb);
                        cut.setUb(ub);
                        cut.setRow(indices.size(), indices.data(), values.data());
                        if (cut.violated(s) > 1.0e-5)
                        {
                            cs.insert(cut);
                            addedWeakCuts = true;
                            cutcntweak++;
                        }
                    }
                }
                else // if (instance->CUV[u][v] != instance->CUV[v][u] || instance->CUV[u][w] != instance->CUV[w][u])
                {
                    auto kk = i;
                    auto ii = j;
                    auto jj = k;
                    if (kk < ii)
                    {
                        swap(ii, kk);
                    }
                    if (kk < jj)
                    {
                        swap(jj, kk);
                    }
                    auto uu = subgraph[ii];
                    auto vv = subgraph[jj];
                    auto ww = subgraph[kk];
                    vector<int> indices;
                    vector<double> values;
                    int lb = 0;
                    int ub = 1;
                    if (instance->pomatrix[uu][vv] == 0)
                    {
                        indices.push_back(var_index[ii][jj]);
                        values.push_back(1.0);
                    }
                    else if (instance->pomatrix[uu][vv] == 1)
                    {
                        lb -= 1;
                        ub -= 1;
                    }
                    if (instance->pomatrix[vv][ww] == 0)
                    {
                        indices.push_back(var_index[jj][kk]);
                        values.push_back(1.0);
                    }
                    else if (instance->pomatrix[vv][ww] == 1)
                    {
                        lb -= 1;
                        ub -= 1;
                    }
                    if (instance->pomatrix[uu][ww] == 0)
                    {
                        indices.push_back(var_index[ii][kk]);
                        values.push_back(-1.0);
                    }
                    else if (instance->pomatrix[uu][ww] == 1)
                    {
                        lb += 1;
                        ub += 1;
                    }
                    if (indices.size() == 0)
                        continue;
                    OsiRowCut cut;
                    cut.setLb(lb);
                    cut.setUb(ub);
                    cut.setRow(indices.size(), indices.data(), values.data());
                    if (cut.violated(s) > 1.0e-5)
                    {
                        cs.insert(cut);
                        addedWeakCuts = true;
                        strongcutcnt++;
                    }
                }
            }
        }
    }
    debug(cutcntweak);
    debug(strongcutcnt);
}

vector<u_int32_t> solve_ilp(OLCMInstance &instance, vector<u_int32_t> &subgraph, vector<u_int32_t> &multiplicities)
{
    CoinModel build;
    debug("Solving ILP witn number of nodes", subgraph.size());

    static vector<u_int32_t> index(instance.nL2);
    vector<vector<u_int32_t>> var_index(subgraph.size(), vector<u_int32_t>(subgraph.size(), 0));
    vector<pair<u_int32_t, u_int32_t>> var_vals;

    for (u_int32_t i = 0; i < subgraph.size(); i++)
        index[subgraph[i]] = i;

    int32_t lb = 0;
    int cnt = 0;
    for (u_int32_t i = 0; i < subgraph.size(); i++)
    {
        for (u_int32_t j = i + 1; j < subgraph.size(); j++)
        {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (instance.pomatrix[u][v] == 0)
            {
                // we do not know the order
                var_index[i][j] = cnt++;
                var_vals.push_back({i, j});
                int32_t mult = ((int32_t)multiplicities[i]) * ((int32_t)multiplicities[j]);
                auto obj = mult * ((double)instance.CUV[u][v]) - ((double)instance.CUV[v][u]);
                string s = to_string(i) + "_" + to_string(j);
                build.addColumn(0, NULL, NULL, 0, 1, mult * (((double)instance.CUV[u][v]) - ((double)instance.CUV[v][u])), s.c_str(), true);
                lb += ((int32_t)instance.CUV[v][u]) * mult;
            }
        }
    }

    CutInfo info;
    info.var_index = &var_index;
    info.subgraph = &subgraph;
    info.index = &index;
    info.instance = &instance;
    CglTriangles triangles(&info);
    triangles.setGlobalCuts(true);

    int indices[] = {0};
    double values[] = {1.0};
    // empty model is buggy
    build.addRow(1, indices, values, 0, 1);

    string collector;

    OsiClpSolverInterface solver1;

#ifdef MYLOCAL
    solver1.setLogLevel(0);
#else
    BufferStdout swapstdout(collector);
    solver1.setLogLevel(0);
#endif

    solver1.loadFromCoinModel(build);

    CbcModel model(solver1);
    OsiBabSolver defaultC;
    defaultC.setSolverType(4);
#ifdef MYLOCAL
    model.setLogLevel(1);
#else
    model.setLogLevel(0);
#endif
    model.setNumberThreads(1);

    debug(model.getNumRows());

    int n = subgraph.size();
    int nchoose3 = (n * (n - 1) * (n - 2)) / 6;
    debug(nchoose3);

    model.addCutGenerator(&triangles, 1, "lazy", true, 1);
    model.cutGenerator(0)->setMustCallAgain(true);
    HeuristicInfo hinfo;
    hinfo.index = &index;
    hinfo.instance = &instance;
    hinfo.subgraph = &subgraph;
    hinfo.var_index = &var_index;
    hinfo.var_vals = &var_vals;
    MyCbcRoundingHeuristic heuristic(model, &hinfo);
    model.addHeuristic(&heuristic, "myrounding", 3);

    OsiClpSolverInterface *osiclp = dynamic_cast<OsiClpSolverInterface *>(model.solver());
    // go faster stripes
    if (osiclp)
    {
        // Don't allow dual stuff
        osiclp->setSpecialOptions(osiclp->specialOptions() | 262144);
    }

    model.initialSolve();
    model.solver()->setAuxiliaryInfo(&defaultC);
    model.passInSolverCharacteristics(&defaultC);
    // Say don't recompute solution d)
    model.setSpecialOptions(4);

    model.branchAndBound(3);
    instance.LB += lb + model.getSolverObjValue();
    auto s = model.getColSolution();
    vector<u_int32_t> ans;
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
            else // if (instance.CUV[u][v] != instance.CUV[v][u])
            {
                auto val = s[var_index[i][j]];
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
    return ans;
}