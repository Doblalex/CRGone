#pragma once
#include "util.hpp"
#include "instance.hpp"

u_int64_t cuv(u_int32_t u, u_int32_t v, OLCMInstance& instance);

void compute_cuv(OLCMInstance& instance);

class TrieNode {
public:
    std::unordered_map<u_int32_t, shared_ptr<TrieNode>> children;
    std::vector<u_int32_t> holding;
};
     
void reduce_twins(OLCMInstance& instance);

void fix_partial_order(OLCMInstance& instance, vector<u_int32_t>& subgraph);
    
void beginendreduction(vector<u_int32_t>& prefix, vector<u_int32_t>& suffix, vector<u_int32_t>& subgraph, OLCMInstance& instance);

void cuvtwins(OLCMInstance& instance, vector<u_int32_t>& subgraph);