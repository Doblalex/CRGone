#include "decomposition.hpp"
#include "partition.hpp"

vector<vector<u_int32_t>> get_refined_by(u_int32_t node, OLCMInstance &instance, vector<u_int32_t> subgraph)
{
    vector<u_int32_t> partitionel;
    for (auto x : subgraph)
    {
        if (x != node)
            partitionel.push_back(x);
    }
    PartitionRefinement mypartition(partitionel, node);
    while (mypartition.dorefine(instance))
    {
    }

    vector<vector<u_int32_t>> ans;
    ans.emplace_back();
    ans[0].push_back(node);
    auto pel = mypartition.start;
    int cnt = 0;
    bool somegreater = false;
    vector<PartitionElement *> elements;
    while (pel)
    {
        elements.push_back(pel);
        vector<u_int32_t> ansrow;
        for (auto x : pel->elements)
        {
            ansrow.push_back(x);
        }

        ans.push_back(ansrow);
        if (pel->elements.size() > 1)
        {
            debug(pel->elements.size());
            somegreater = true;
        }

        cnt += pel->elements.size();
        pel = pel->next;
    }
    // for (auto module:ans) {
    //     unordered_set<u_int32_t> module_set(module.begin(), module.end());
    //     for (auto u: module) {
    //         for (auto v: module) {
    //             for (auto x: subgraph) {
    //                 if (module_set.find(x) == module_set.end()) {
    //                     assert(instance.CUV[u][x] == instance.CUV[v][x] && instance.CUV[x][u] == instance.CUV[x][v]);
    //                 }
    //             }
    //         }
    //     }
    // }
    for (u_int32_t i = 0; i < ans.size(); i++)
    {
        auto module1 = ans[i];
        for (u_int32_t j = i + 1; j < ans.size(); j++)
        {
            auto module2 = ans[j];
            bool forward = false;
            bool backward = false;
            for (auto u : module1)
            {
                for (auto v : module2)
                {
                    assert(u != v);
                    if (instance.pomatrix[u][v] == 1) {
                        forward = true;
                    }
                        
                    if (instance.pomatrix[v][u] == 1) {
                        backward = true;
                    }
                        
                }
            }
            assert(!backward || !forward);
            if (backward || forward)
            {
                for (auto u : module1)
                {
                    for (auto v : module2)
                    {
                        if (forward && instance.pomatrix[u][v] != 1)
                        {
                            instance.pomatrix[u][v] = 1;
                            instance.pomatrix[v][u] = -1;
                            instance.LB += instance.CUV[u][v];
                        }
                        if (backward && instance.pomatrix[v][u] != 1) {
                            instance.pomatrix[v][u] = 1;
                            instance.pomatrix[u][v] = -1;
                            instance.LB += instance.CUV[v][u];
                        }
                    }
                }
            }
        }
    }
    assert(cnt + 1 == subgraph.size());
    return ans;
}