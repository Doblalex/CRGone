#pragma once

#include "util.hpp"

class OLCMInstance
{
public:
    u_int32_t nL1;
    u_int32_t nL2;
    u_int32_t m;
    vector<vector<u_int32_t>> al;
    vector<vector<u_int64_t>> CUV;
    vector<vector<int8_t>> pomatrix;
    vector<vector<u_int32_t>> twincorrespondents;
    vector<u_int32_t> working_nodes;
    u_int64_t LB;

    OLCMInstance() {}
    ~OLCMInstance() = default;

    OLCMInstance(u_int32_t nL1, u_int32_t nL2, u_int64_t m, vector<vector<u_int32_t>> al) : nL1(nL1), nL2(nL2), m(m), al(al)
    {
        LB = 0;
    }

    OLCMInstance(OLCMInstance& t)
    {
        nL1 = t.nL1;
        nL2 = t.nL2;
        m = t.m;
        al = t.al;
        CUV = t.CUV;
        pomatrix = t.pomatrix;
        twincorrespondents = t.twincorrespondents;
        working_nodes = t.working_nodes;
        LB = t.LB;
    }
};