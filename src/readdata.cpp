#include "readdata.hpp"

OLCMInstance read_instance(bool param) {
    u_int32_t nL1, nL2;
    u_int32_t m;
    string s;
    string line;
    getline(cin, line);
    istringstream iss(line);
    iss>>s;
    iss>>s;
    iss >> nL1 >> nL2 >> m;
    int abc;
    if (iss>>abc) {
        param = true;
    }

    if (param) {
        for (int i = 0; i < nL1 + nL2; i++) cin>>s;
        getline(cin, line);
    }    

    std::vector<std::vector<u_int32_t>> adj(nL2);

    while (m--) {
        std::getline(std::cin, line);
        if (line.find("c ") != std::string::npos || line.empty())
            continue;
        u_int32_t u1, u2;
        std::istringstream iss(line);
        iss >> u1 >> u2;

        adj[u2 - nL1 - 1].push_back(u1 - 1);
    }

    for (int u = 0; u < nL2; ++u) {
        std::sort(adj[u].begin(), adj[u].end());
    }
    return OLCMInstance(nL1, nL2, m, adj);
}