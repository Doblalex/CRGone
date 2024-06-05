// A C++ program to find biconnected components in a given undirected graph
#include "scc.hpp"

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v, int w) { adj[v].push_back(w); }

// A recursive function that finds and prints strongly
// connected components using DFS traversal u --> The vertex
// to be visited next disc[] --> Stores discovery times of
// visited vertices low[] -- >> earliest visited vertex (the
// vertex with minimum
//             discovery time) that can be reached from
//             subtree rooted with current vertex
// *st -- >> To store all the connected ancestors (could be
// part
//         of SCC)
// stackMember[] --> bit/index array for faster check
// whether
//                 a node is in stack
void Graph::SCCUtil(int u, int disc[], int low[],
                           stack<int> *st, bool stackMember[], vector<vector<int>>& sccs)
{
    // A static variable is used for simplicity, we can
    // avoid use of static variable by passing a pointer.
    static int time = 0;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    st->push(u);
    stackMember[u] = true;

    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i; // v is current adjacent of 'u'

        // If v is not visited yet, then recur for it
        if (disc[v] == -1)
        {
            SCCUtil(v, disc, low, st, stackMember, sccs);

            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            // Case 1 (per above discussion on Disc and Low
            // value)
            low[u] = min(low[u], low[v]);
        }

        // Update low value of 'u' only of 'v' is still in
        // stack (i.e. it's a back edge, not cross edge).
        // Case 2 (per above discussion on Disc and Low
        // value)
        else if (stackMember[v] == true)
            low[u] = min(low[u], disc[v]);
    }

    // head node found, pop the stack and print an SCC
    int w = 0; // To store stack extracted vertices
    vector<int> ans;
    if (low[u] == disc[u])
    {
        while (st->top() != u)
        {
            w = (int)st->top();
            ans.push_back(w);
            stackMember[w] = false;
            st->pop();
        }
        w = (int)st->top();
        ans.push_back(w),
            stackMember[w] = false;
        st->pop();
    }
    if (ans.size() > 0) {
        sccs.push_back(ans);
    }
}

// The function to do DFS traversal. It uses SCCUtil()
vector<vector<int>> Graph::SCC()
{
    int *disc = new int[V];
    int *low = new int[V];
    bool *stackMember = new bool[V];
    stack<int> *st = new stack<int>();

    // Initialize disc and low, and stackMember arrays
    for (int i = 0; i < V; i++)
    {
        disc[i] = NIL;
        low[i] = NIL;
        stackMember[i] = false;
    }

    // Call the recursive helper function to find strongly
    // connected components in DFS tree with vertex 'i'
    vector<vector<int>> ans;
    for (int i = 0; i < V; i++)
    {
        if (disc[i] == NIL)
        {
            SCCUtil(i, disc, low, st, stackMember, ans);
        }
    }

    return ans;
}