#include "preprocess.hpp"

u_int64_t cuv(u_int32_t u, u_int32_t v, OLCMInstance& instance) {
    u_int32_t iv = 0;
    u_int64_t ans = 0;
    if (instance.al[u].size() == 0 || instance.al[v].size() == 0) {
        return 0;
    }
    if (instance.al[u].back() <= instance.al[v].front()) {
        return 0;
    } else if (instance.al[u].front() > instance.al[v].back()) {
        return instance.al[u].size() * instance.al[v].size();
    }
    for (int x : instance.al[u]) {
        while (iv < instance.al[v].size() && instance.al[v][iv] < x) {
            iv++;
        }
        ans += iv;
    }
    return ans;
}

void compute_cuv(OLCMInstance& instance) {
    for (int u : instance.working_nodes) {
        for (int v : instance.working_nodes) {
            instance.CUV[u][v] = cuv(u, v, instance) * instance.twincorrespondents[u].size() * instance.twincorrespondents[v].size();
        }
    }
    for (u_int32_t i = 0; i < instance.working_nodes.size(); i++) {
        for (u_int32_t j = i+1; j < instance.working_nodes.size(); j++) {
            auto u = instance.working_nodes[i];
            auto v = instance.working_nodes[j];
            auto mn = min(instance.CUV[u][v], instance.CUV[v][u]);
            instance.LB += mn;
            instance.CUV[u][v] -= mn;
            instance.CUV[v][u] -= mn;
        }
    }
}

void reduce_twins(OLCMInstance& instance) {
    std::vector<u_int32_t> nodesstay;
    instance.twincorrespondents.resize(instance.nL2);
    
    auto root = std::make_shared<TrieNode>();
    std::vector<std::shared_ptr<TrieNode>> nodes;
    nodes.push_back(root);

    for (u_int32_t u = 0; u < instance.nL2; ++u) {
        auto node = root;
        for (u_int32_t v : instance.al[u]) {
            if (node->children.find(v) == node->children.end()) {
                node->children[v] = std::make_shared<TrieNode>();
                nodes.push_back(node->children[v]);
            }
            node = node->children[v];
        }
        node->holding.push_back(u);
    }

    for (const auto& node : nodes) {
        if (!node->holding.empty()) {
            nodesstay.push_back(node->holding[0]);
            instance.twincorrespondents[node->holding[0]] = node->holding;
            instance.LB += (node->holding.size() * (node->holding.size() - 1) * cuv(node->holding[0], node->holding[0], instance)) / 2;
        }
    }

    std::sort(nodesstay.begin(), nodesstay.end());
    instance.working_nodes = nodesstay;
}

void fix_partial_order(OLCMInstance& instance, vector<u_int32_t>& subgraph) {
    if (subgraph.size() == 0) return;
    std::unordered_set<std::pair<u_int32_t, u_int32_t>, pair_hash> ordered;
    debug("doing median reduction for #nodes", subgraph.size());
    for (u_int32_t i = 0; i < subgraph.size(); ++i) {
        for (u_int32_t j = i + 1; j < subgraph.size(); ++j) {
            u_int32_t u = subgraph[i];
            u_int32_t v = subgraph[j];
            if (instance.al[u][instance.al[u].size()-1] <= instance.al[v][0]){
                instance.pomatrix[u][v] = 1;
                instance.pomatrix[v][u] = -1;
                continue;
            } else if (instance.al[v][instance.al[v].size()-1] <= instance.al[u][0]) {
                instance.pomatrix[v][u] = 1;
                instance.pomatrix[u][v] = -1;
                continue;
            }

            int swapcnt = 1;
            while (swapcnt--) {
                
                if (instance.al[u].size() == 2 && instance.al[v].size() == 2) {
                    if (instance.al[u][0] == instance.al[v][0] && instance.al[u][1] > instance.al[v][1]) {
                        debug("RR3 case 1 applies");
                        instance.pomatrix[v][u] = 1;
                        instance.pomatrix[u][v] = -1;
                        instance.LB += instance.CUV[v][u];
                    }
                    else if (instance.al[u][1] == instance.al[v][1] && instance.al[u][0] > instance.al[v][0]) {
                        debug("RR3 case 2 applies");
                        instance.pomatrix[v][u] = 1;
                        instance.pomatrix[u][v] = -1;
                        instance.LB += instance.CUV[v][u];
                    }
                }
                swap(u,v);
            }
        }
    }    
    // for (u_int32_t i = 0; i < subgraph.size(); i++)
    // {
    //     for (u_int32_t j = 0; j < subgraph.size(); j++)
    //     {
    //         for(u_int32_t k = 0; k < subgraph.size(); k++)
    //         {
    //             auto u = subgraph[i];
    //             auto v = subgraph[j];
    //             auto w = subgraph[k];
    //             if (instance.pomatrix[u][v] == 1 && instance.pomatrix[v][w] == 1) {
    //                 assert(instance.pomatrix[u][w] == 1);
    //             }
    //         }
    //     }
    // }
}

void beginendreduction(vector<u_int32_t>& prefix, vector<u_int32_t>& suffix, vector<u_int32_t>& subgraph, OLCMInstance& instance) {
    static vector<u_int32_t> betterbeforecnt(instance.nL2, 0);
    static vector<u_int32_t> betteraftercnt(instance.nL2, 0);

    for (auto u: subgraph) {
        betterbeforecnt[u] = 0;
        betteraftercnt[u] = 0;
    }

    for (u_int32_t i = 0; i < subgraph.size(); i++) {
        for (u_int32_t j = i+1; j < subgraph.size(); j++) {
            auto u = subgraph[i];
            auto v = subgraph[j];
            if (instance.CUV[u][v] == instance.CUV[v][u]) {
                betterbeforecnt[u]++;
                betterbeforecnt[v]++;
                betteraftercnt[u]++;
                betteraftercnt[v]++;
            }
            else if(instance.CUV[u][v] < instance.CUV[v][u]) {
                betterbeforecnt[u]++;
                betteraftercnt[v]++;
            }
            else {
                betteraftercnt[u]++;
                betterbeforecnt[v]++;
            }
        }
    }

    unordered_set<u_int32_t> nodes(subgraph.begin(), subgraph.end());
    bool progress = true;
    while (progress) {
        progress = false;
        int nodebegin = -1;
        int nodeend = -1;
        for (auto u: nodes) {
            if (betterbeforecnt[u] == nodes.size() - 1) {
                nodebegin = u;
                break;
            }
            if (betteraftercnt[u] == nodes.size() - 1) {
                nodeend = u;
                break;
            }
        }
        if (nodeend != -1 || nodebegin != -1) {
            progress = true;
            u_int32_t u = 0;
            if (nodebegin != -1) {
                u = nodebegin;
                for (auto v: nodes) {
                    if (u != v && instance.pomatrix[u][v] != 1) {
                        instance.LB += instance.CUV[u][v];
                    }
                        
                }
                prefix.push_back(nodebegin);                
            }
            else {
                u = nodeend;
                for (auto v: nodes) {
                    if (u != v && instance.pomatrix[v][u] != 1) {
                        instance.LB += instance.CUV[v][u];
                    }
                        
                }
                suffix.push_back(nodeend);                
            }
            nodes.erase(u);
            for (auto v: nodes) {
                if (instance.CUV[u][v] == instance.CUV[v][u]) {
                    betterbeforecnt[v]--;
                    betteraftercnt[v]--;
                }
                else if(instance.CUV[u][v] < instance.CUV[v][u]) {
                    betteraftercnt[v]--;
                }
                else {
                    betterbeforecnt[v]--;
                }
            }
        }
    }
    subgraph = vector<u_int32_t>(nodes.begin(), nodes.end());
    std::sort(subgraph.begin(), subgraph.end());
    std::reverse(suffix.begin(), suffix.end());
}

void cuvtwins(OLCMInstance& instance, vector<u_int32_t>& subgraph) {
    int cnt = 0;
    for (auto u: subgraph) {
        for (auto v: subgraph) {
            if (u<v) {
                bool istwins = true;
                for (auto w: subgraph) {
                    if (u != w && v != w && (instance.CUV[u][w]-instance.CUV[w][u]) != (instance.CUV[v][w]-instance.CUV[w][v])) {
                        istwins = false;
                        break;
                    }
                }
                if (istwins) {
                    cnt++;
                }
            }
        }
    }
    cout<<"twincnt "<<cnt<<endl;
}