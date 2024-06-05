#include "mycbcheuristic.hpp"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>

#include "CoinHelperFunctions.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinTime.hpp"

// Default Constructor
MyCbcRoundingHeuristic::MyCbcRoundingHeuristic()
    : CbcHeuristic()
{
}

// Constructor with model - assumed before cuts
MyCbcRoundingHeuristic::MyCbcRoundingHeuristic(CbcModel &model, HeuristicInfo *info)
    : CbcHeuristic(model)
{
    this->myInfo = info;
    model_ = &model;
    incoming = vector<u_int32_t>(myInfo->subgraph->size(), 0);
    g = vector<vector<u_int32_t>>(myInfo->subgraph->size());
    for (u_int32_t i = 0; i < myInfo->subgraph->size(); i++) {
        for (u_int32_t j = i+1; j < myInfo->subgraph->size(); j++) {
            auto val = myInfo->instance->pomatrix[(*myInfo->subgraph)[i]][(*myInfo->subgraph)[j]];
            if (val == 1) {
                g[i].push_back(j);
                incoming[j]++;
            }
            else if(val == -1) {
                incoming[i]++;
                g[j].push_back(i);
            }
        }
    }
    setWhen(3);
}

// Destructor
MyCbcRoundingHeuristic::~MyCbcRoundingHeuristic()
{
}

// Clone
CbcHeuristic *
MyCbcRoundingHeuristic::clone() const
{
    return new MyCbcRoundingHeuristic(*this);
}

// Copy constructor
MyCbcRoundingHeuristic::MyCbcRoundingHeuristic(const MyCbcRoundingHeuristic &rhs)
    : CbcHeuristic(rhs)
{
    this->myInfo = rhs.myInfo;
    this->g = rhs.g;
    this->incoming = rhs.incoming;
}

// Assignment operator
MyCbcRoundingHeuristic &
MyCbcRoundingHeuristic::operator=(const MyCbcRoundingHeuristic &rhs)
{
    if (this != &rhs)
    {
        CbcHeuristic::operator=(rhs);
    }
    return *this;
}

// Resets stuff if model changes
void MyCbcRoundingHeuristic::resetModel(CbcModel *model)
{
    CbcHeuristic::resetModel(model);
}

bool MyCbcRoundingHeuristic::shouldHeurRun(int whereFrom) {
    debug("in should heuristic run");
    return true;
}

/*
  Randomized Rounding Heuristic
  Returns 1 if solution, 0 if not
*/
int MyCbcRoundingHeuristic::solution(double &solutionValue,
                                     double *betterSolution)
{
    auto s = model_->getColSolution();
    vector<u_int32_t> t_incoming(incoming.begin(), incoming.end());
    vector<double> sortby(myInfo->subgraph->size(), 0.0);
    double solfixedvars = 0;
    for (u_int32_t i = 0; i < myInfo->subgraph->size(); i++) {
        sortby[i] += t_incoming[i];
    }
    for (u_int32_t i = 0; i < myInfo->var_vals->size(); i++) {
        auto p = (*myInfo->var_vals)[i];
        sortby[p.first] += (1.0-s[i]);
        sortby[p.second] += s[i];
    }
    set<pair<double,u_int32_t>> free;
    for (u_int32_t i = 0; i < myInfo->subgraph->size(); i++) {
        if (t_incoming[i] == 0) {
            free.insert({sortby[i], i});
        }
    }
    vector<u_int32_t> order;
    unordered_map<u_int32_t, u_int32_t> order_index;
    int cnt = 0;
    while (free.size()) {
        auto p = *free.begin();
        auto i = p.second;
        free.erase(p);
        order.push_back(i);
        order_index[i] = cnt;
        cnt++;
        for (auto j: g[i]) {
            if (--t_incoming[j] == 0) {
                free.insert({sortby[j], j});
            }
        }
    }
    
    double newsolutionvalue = 0;
    for (u_int32_t x = 0; x < myInfo->var_vals->size(); x++) {
        auto p = (*myInfo->var_vals)[x];
        auto i = p.first;
        auto j = p.second;
        auto u = (*myInfo->subgraph)[i];
        assert(i != j);
        if (order_index[i]<order_index[j]) {
            betterSolution[x] = 1.0;
            newsolutionvalue+=model_->getObjCoefficients()[x];
        }
        else {
            betterSolution[x] = 0.0;
        }
    }
    debug(newsolutionvalue);
    if (newsolutionvalue < solutionValue) {
        debug("new best solution", newsolutionvalue);
        // model_->setBestSolution(betterSolution, myInfo->var_vals->size(), COIN_DBL_MAX, true);
        solutionValue = newsolutionvalue;
        return 1;
    }

    return 0;
}

int MyCbcRoundingHeuristic::solution2(double &solutionValue,
                                     double *betterSolution)
{
    debug("in rounding heuristic2");
    return 0;
}
// update model
void MyCbcRoundingHeuristic::setModel(CbcModel *model)
{
    CbcHeuristic::setModel(model);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
