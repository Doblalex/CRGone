#pragma once
#include "util.hpp"
#include "instance.hpp"
#include "CbcHeuristic.hpp"

typedef struct
{
    vector<vector<u_int32_t>> *var_index;
    vector<u_int32_t> *subgraph;
    OLCMInstance *instance;
    vector<u_int32_t> *index;
    vector<pair<u_int32_t, u_int32_t>>* var_vals;
} HeuristicInfo;

class MyCbcRoundingHeuristic : public CbcHeuristic
{
public:
    HeuristicInfo* myInfo;
    vector<u_int32_t> incoming;
    vector<vector<u_int32_t>> g;
    // Default Constructor
    MyCbcRoundingHeuristic();

    /* Constructor with model - assumed before cuts
         Initial version does not do Lps
      */
    MyCbcRoundingHeuristic(CbcModel &model, HeuristicInfo *info);

    // Copy constructor
    MyCbcRoundingHeuristic(const MyCbcRoundingHeuristic &);

    // Destructor
    ~MyCbcRoundingHeuristic();

    /// Clone
    virtual CbcHeuristic *clone() const;

    /// Assignment operator
    MyCbcRoundingHeuristic &operator=(const MyCbcRoundingHeuristic &rhs);

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel *model);

    /// update model (This is needed if cliques update matrix etc)
    virtual void setModel(CbcModel *model);

    using CbcHeuristic::solution;
    /** returns 0 if no solution, 1 if valid solution.
          Sets solution values if good, sets objective value (only if good)
          needs comments
      */
    virtual int solution(double &objectiveValue,
                         double *newSolution);

    virtual int solution2(double &objectiveValue,
                         double *newSolution);

    virtual bool shouldHeurRun(int whereFrom);

protected:
};