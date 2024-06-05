#pragma once

#include "util.hpp"
#include "instance.hpp"

bool separate_cut(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<u_int32_t>& left, vector<u_int32_t>& right);

// bool separate_bcc(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<vector<u_int32_t>>& partition);

void separate_scc(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<vector<u_int32_t>>& partition);

vector<vector<u_int32_t>> separate_adjacency(OLCMInstance& instance, vector<u_int32_t>& subgraph);

vector<vector<u_int32_t>> separate_partial_order(vector<unordered_set<u_int32_t>>& po, vector<unordered_set<u_int32_t>>& porev, vector<u_int32_t>& subgraph);