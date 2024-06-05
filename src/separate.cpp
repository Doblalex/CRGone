#include "separate.hpp"
#include "scc.hpp"

// bool separate_bcc(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<vector<u_int32_t>>& partition) {
//     Graph g(subgraph.size());
//     for (u_int32_t i = 0; i < subgraph.size(); i++) {
//         for (u_int32_t j = i+1; j < subgraph.size(); j++) {
//             auto u = subgraph[i];
//             auto v = subgraph[j];
//             if (instance.CUV[u][v] != instance.CUV[v][u]) {
//                 g.addEdge(i, j);
//                 g.addEdge(j, i);
//             }
//         }
//     }
//     g.BCC();
//     return false;
// }

void separate_scc(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<vector<u_int32_t>>& partition) {
    Graph g(subgraph.size());
    for (u_int32_t i = 0; i < subgraph.size(); i++) {
        for (u_int32_t j = i+1; j < subgraph.size(); j++) {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (instance.CUV[u][v] > instance.CUV[v][u]) {
                g.addEdge(j,i);
            }
            else if (instance.CUV[v][u] > instance.CUV[u][v]) {
                g.addEdge(i,j);
            }
        }
    }

    auto sccs = g.SCC();
    for (auto scc: sccs)debug(subgraph.size(), scc.size());
    // need to sort sccs
    vector<int> inscc(subgraph.size());
    for (int i = 0; i < sccs.size(); i++) {
        for (auto x: sccs[i]) {
            inscc[x] = i;
        }
    }

    vector<unordered_set<int>> sccgr(sccs.size());
    vector<unordered_set<int>> sccgrrev(sccs.size());

    for (u_int32_t i = 0; i < subgraph.size(); i++) {
        for (u_int32_t j = i+1; j < subgraph.size(); j++) {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (inscc[i] == inscc[j]) continue;
            if (instance.CUV[u][v] > instance.CUV[v][u]) {
                sccgr[inscc[j]].insert(inscc[i]);
                sccgrrev[inscc[i]].insert(inscc[j]);
            }
            else if (instance.CUV[v][u] > instance.CUV[u][v]) {
                sccgr[inscc[i]].insert(inscc[j]);
                sccgrrev[inscc[j]].insert(inscc[i]);
            }
        }
    }

    vector<int> incoming(sccgr.size(), 0);
    vector<int> s;
    for (int i = 0; i < sccs.size(); i++) {
        incoming[i] = sccgrrev[i].size();
        if (incoming[i] == 0) s.push_back(i);
    }
    partition.clear();
    while (s.size() > 0) {
        auto i = *s.rbegin();
        s.pop_back();
        partition.emplace_back();
        for (auto j: sccgr[i]) {
            incoming[j]--;
            if (incoming[j] == 0) {
                s.push_back(j);
            }
        }
        for (auto x: sccs[i]) {
            partition.rbegin()->push_back(subgraph[x]);
        }
    }
}

unordered_set<u_int32_t> reachable(vector<vector<u_int32_t>>& g, u_int32_t start) {
    unordered_set<u_int32_t> ans;
    ans.insert(start);
    vector<u_int32_t> q;
    q.push_back(start);
    while (q.size()) {
        auto u = *q.rbegin();
        q.pop_back();
        for (auto v: g[u]) {
            if (ans.find(v) == ans.end()) {
                ans.insert(v);
                q.push_back(v);
            }
        }
    }
    return ans;
}

bool separate_cut(OLCMInstance& instance, vector<u_int32_t>& subgraph, vector<u_int32_t>& left, vector<u_int32_t>& right) {
    vector<vector<u_int32_t>> g(subgraph.size());
    vector<vector<u_int32_t>> g_rev(subgraph.size());
    debug(subgraph.size());
    for (u_int32_t i = 0; i < subgraph.size(); i++) {
        for (u_int32_t j = i+1; j < subgraph.size(); j++) {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (instance.CUV[u][v] > instance.CUV[v][u]) {
                g_rev[j].push_back(i);
                g[i].push_back(j);
            }
            else if (instance.CUV[v][u] > instance.CUV[u][v]) {
                g[j].push_back(i);
                g_rev[i].push_back(j);
            }
        }
    }

    // for (u_int32_t i = 0; i < subgraph.size(); i++) {
    //         debug(i, g[i].size(), g_rev[i].size());
    // }
    // exit(0);
    auto a = reachable(g, 0);
    if (a.size() != subgraph.size()) {
        debug("separation", a.size(), subgraph.size());
        for (auto x: a) {
            left.push_back(subgraph[x]);
        }
        for (u_int32_t i = 0; i < subgraph.size(); i++) {
            if (a.find(i) == a.end()) {
                right.push_back(subgraph[i]);
            }
        }
        sort(left.begin(), left.end());
        for (auto u: left) {
            for (auto v: right) {
                instance.LB += instance.CUV[u][v];
            }
        }
        return true;
    }
    a = reachable(g_rev, 0);
    if (a.size() != subgraph.size()) {
        debug("separation", a.size(), subgraph.size());
        for (auto x: a) {
            right.push_back(subgraph[x]);
        }
        for (u_int32_t i = 0; i < subgraph.size(); i++) {
            if (a.find(i) == a.end()) {
                left.push_back(subgraph[i]);
            }
        }
        for (auto u: left) {
            for (auto v: right) {
                instance.LB += instance.CUV[u][v];
            }
        }
        sort(right.begin(), right.end());
        return true;
    }
    return false;
}

vector<vector<u_int32_t>> separate_adjacency(OLCMInstance& instance, vector<u_int32_t>& subgraph) {
    vector<tuple<u_int32_t, int, u_int32_t>> eventqueue;
    vector<vector<u_int32_t>> subgraphs;
    for (auto u: subgraph) {        
        if (instance.al[u].size() == 1) {
            eventqueue.push_back({instance.al[u][instance.al[u].size()-1]*2+1, 0, u});
            eventqueue.push_back({instance.al[u][0]*2, 1, u});
        }
        else {
            eventqueue.push_back({instance.al[u][instance.al[u].size()-1]*2, 0, u});
            eventqueue.push_back({instance.al[u][0]*2, 1, u});
        }
        
    }
    sort(eventqueue.begin(), eventqueue.end());
    unordered_set<int> active;
    vector<u_int32_t> cursubgraph;
    for (auto event: eventqueue) {
        auto u = get<2>(event);
        if (get<1>(event) == 0) {
            active.erase(u);
            if (active.size() == 0) {
                sort(cursubgraph.begin(), cursubgraph.end());
                subgraphs.push_back(cursubgraph);
                cursubgraph = vector<u_int32_t>();
            }
        } else {
            active.insert(u);
            cursubgraph.push_back(u);
        }
    }
    return subgraphs;
}

vector<vector<u_int32_t>> separate_partial_order(vector<unordered_set<u_int32_t>>& po, vector<unordered_set<u_int32_t>>& porev, vector<u_int32_t>& subgraph) {
    std::vector<u_int32_t> nodesstay = subgraph;
    std::vector<std::vector<u_int32_t>> subgraphs;

    while (!nodesstay.empty()) {
        std::sort(nodesstay.begin(), nodesstay.end(), [&po](int u1, int u2) {
            return po[u1].size() > po[u2].size();
        });

        std::unordered_map<u_int32_t, u_int32_t> index;
        for (size_t i = 0; i < nodesstay.size(); ++i)
            index[nodesstay[i]] = i;

        std::unordered_set<u_int32_t> subgraph;
        std::vector<u_int32_t> incoming_cnt(nodesstay.size(), 0);
        u_int32_t outdeg = 0;
        u_int32_t i = 0;
        while (i < nodesstay.size()) {
            u_int32_t u = nodesstay[i];
            outdeg -= incoming_cnt[i];
            subgraph.insert(u);
            for (u_int32_t v : po[nodesstay[i]]) {
                if (index[v] > index[u]) {
                    outdeg++;
                    incoming_cnt[index[v]]++;
                }
            }
            u_int32_t rest = nodesstay.size() - subgraph.size();
            if (outdeg == static_cast<u_int32_t>(rest) * static_cast<u_int32_t>(subgraph.size())) {
                
                subgraphs.push_back(std::vector<u_int32_t>(subgraph.begin(), subgraph.end()));
                std::unordered_set<u_int32_t> restgraph(nodesstay.begin(), nodesstay.end());
                for (auto u: subgraph) restgraph.erase(u);
                for (u_int32_t u : subgraph) {
                    for (auto v: vector<u_int32_t>(po[u].begin(), po[u].end())) {
                        if (restgraph.find(v) != restgraph.end()) {
                            po[u].erase(v);
                            porev[v].erase(u);
                        }
                    }
                    for (auto v: vector<u_int32_t>(porev[u].begin(), porev[u].end())) {
                        if (restgraph.find(v) != restgraph.end()) {
                            porev[u].erase(v);
                            po[v].erase(u);
                        }
                    }
                }
                nodesstay.erase(std::remove_if(nodesstay.begin(), nodesstay.end(), [&subgraph](int u) { return subgraph.count(u) > 0; }), nodesstay.end());
                break;
            }
            i++;
        }
    }
    return subgraphs;
}