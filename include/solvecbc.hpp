#pragma once

#include "util.hpp"
#include "instance.hpp"

vector<u_int32_t> solve_ilp (OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<u_int32_t>& multiplicities);