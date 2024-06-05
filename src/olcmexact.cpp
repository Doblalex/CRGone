#include "util.hpp"
#include "instance.hpp"
#include "preprocess.hpp"
#include "readdata.hpp"
#include "separate.hpp"
#include "solvescip.hpp"
#include "decomposition.hpp"
#include "cxxopts.hpp"

vector<u_int32_t> recursive_solve(OLCMInstance &instance, vector<u_int32_t> subgraph, bool computed_po)
{
    vector<vector<u_int32_t>> subgraphs;
    vector<u_int32_t> ans;
    if (subgraph.size() == 0) return ans;
    if (subgraph.size() == 1) {
        for (auto v : instance.twincorrespondents[subgraph[0]])
        {
            ans.push_back(v);
        }
        return ans;
    }
    subgraphs = separate_adjacency(instance, subgraph);
    if (subgraphs.size() > 1)
    {
        for (auto subgraphinner : subgraphs)
        {
            for (auto u : recursive_solve(instance, subgraphinner, computed_po))
            {
                ans.push_back(u);
            }
        }
        return ans;
    }

    vector<vector<u_int32_t>> partition;
    separate_scc(instance, subgraph, partition);
    if (partition.size() > 1)
    {
        for (auto subgraphinner : partition)
        {
            for (auto u : recursive_solve(instance, subgraphinner, computed_po))
            {
                ans.push_back(u);
            }
        }
        return ans;
    }
    vector<u_int32_t> prefix;
    vector<u_int32_t> suffix;
    beginendreduction(prefix, suffix, subgraph, instance);
    if (prefix.size() > 0 || suffix.size() > 0)
    {
        for (auto u : prefix)
        {
            for (auto v : instance.twincorrespondents[u])
            {
                ans.push_back(v);
            }
        }
        for (auto v : recursive_solve(instance, subgraph, computed_po))
        {
            ans.push_back(v);
        }
        for (auto u : suffix)
        {
            for (auto v : instance.twincorrespondents[u])
            {
                ans.push_back(v);
            }
        }
        return ans;
    }
    // vector<vector<u_int32_t>> partition;
    // separate_bcc(instance, subgraph, partition);
    if (!computed_po)
    {
        fix_partial_order(instance, subgraph);
        computed_po = true;
    }
    // subgraphs = separate_partial_order(po, porev, subgraph);
    // if (subgraphs.size() > 1)
    // {
    //     for (auto subgraphinner : subgraphs)
    //     {
    //         for (auto v : recursive_solve(instance, subgraphinner, po, porev, computed_po))
    //         {
    //             ans.push_back(v);
    //         }
    //     }
    //     return ans;
    // }

    if (subgraph.size())
    {
        
        vector<u_int32_t> nodesilp;
        vector<u_int32_t> multiplicities;
        unordered_map<u_int32_t, vector<u_int32_t>> orders;
        // for (auto u: subgraph) {
        //     nodesilp.push_back(u);
        //     multiplicities.push_back(1);
        //     vector<u_int32_t> o;
        //     o.push_back(u);
        //     orders[u] = recursive_solve(instance, o, computed_po);
        // }
        auto x = get_refined_by(subgraph[0], instance, subgraph);
        vector<pair<u_int32_t, u_int32_t>> ilpnodes_mult;
        for (auto y : x)
        {
            ilpnodes_mult.push_back({y[0], y.size()});
            orders[y[0]] = recursive_solve(instance, y, computed_po);
        }
        sort(ilpnodes_mult.begin(), ilpnodes_mult.end());
        for (auto p: ilpnodes_mult) {
            nodesilp.push_back(p.first);
            multiplicities.push_back(p.second);
        }
        for (auto u : solve_ilp(instance, nodesilp, multiplicities))
        {
            for (auto v : orders[u])
            {
                ans.push_back(v);
            }
        }
    }

    return ans;
}

int main(int argc, char** argv)
{
    // call("gr21.dist");
    OLCMInstance instance = read_instance();
    reduce_twins(instance);

    vector<u_int32_t> prefix;
    for (auto u : instance.working_nodes)
    {
        if (instance.al[u].size() == 0)
        {
            for (auto v : instance.twincorrespondents[u])
            {
                prefix.push_back(v);
            }
            instance.working_nodes.erase(std::remove(instance.working_nodes.begin(), instance.working_nodes.end(), u), instance.working_nodes.end());
        }
    }
    instance.CUV = vector<vector<u_int64_t>>(instance.nL2, std::vector<u_int64_t>(instance.nL2, 0));
    instance.pomatrix = vector<vector<int8_t>>(instance.nL2, std::vector<int8_t>(instance.nL2, 0));
    compute_cuv(instance);

    auto order = recursive_solve(instance, instance.working_nodes, false);
    debug(instance.LB);
    #ifdef MYLOCAL
    int32_t verifyobj = 0;
    for (u_int32_t i = 0; i < order.size(); i++) {
        for (u_int32_t j = i+1; j < order.size(); j++) {
            verifyobj += cuv(order[i], order[j], instance);
        }
    }
    debug(verifyobj);
    assert(verifyobj == instance.LB);
    #endif
    #ifndef MYLOCAL
    for (auto u : prefix)
    {
        cout << u + instance.nL1 + 1 << endl;
    }
    for (auto u : order)
    {
        cout << u + instance.nL1 + 1 << endl;
    }
    #endif
}