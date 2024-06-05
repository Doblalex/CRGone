#include "partition.hpp"

PartitionRefinement::PartitionRefinement(vector<u_int32_t> V, u_int32_t v)
{
    int n = V.size();
    this->start = new PartitionElement();
    this->end = this->start;
    for (int i = 0; i < n; i++)
    {
        this->start->elements.push_back(V[i]);
        this->start->needrefineby.push_back(v);
    }
    nonemptyrefineby.push_back(this->start);
}

bool PartitionRefinement::dorefine(OLCMInstance &instance)
{
    if (nonemptyrefineby.size() == 0)
        return false;
    auto x = nonemptyrefineby.front();
    nonemptyrefineby.pop_front();
    if (x->elements.size() == 1) {
        return true;
    }
        
    

    u_int32_t u = x->needrefineby[0];

    unordered_map<int, vector<u_int32_t>> bycol;
    for (auto v : x->elements)
    {
        int cuvdiff = ((int)instance.CUV[u][v]) - ((int)instance.CUV[v][u]);
        bycol[cuvdiff].push_back(v);
    }

    vector<PartitionElement *> newelements;
    for (auto &it : bycol)
    {
        PartitionElement *newel = new PartitionElement();
        for (auto v : it.second)
        {
            newel->elements.push_back(v);
        }
        for (auto uu : x->needrefineby)
        {
            if (uu != u)
                newel->needrefineby.push_back(uu);
        }
        newel->last = this->end;
        this->end->next = newel;
        this->end = newel;
        newelements.push_back(newel);
    }
    for (u_int32_t i = 0; i < newelements.size(); i++)
    {
        for (u_int32_t j = i + 1; j < newelements.size(); j++)
        {
            auto el1 = newelements[i];
            auto el2 = newelements[j];
            for (auto uu : el1->elements)
            {
                el2->needrefineby.push_back(uu);
            }
            for (auto uu : el2->elements)
            {
                el1->needrefineby.push_back(uu);
            }
        }
    }

    for (auto el: newelements) {
        if (el->needrefineby.size() > 0) nonemptyrefineby.push_back(el);
    }

    if (x->next == NULL)
    {
        x->last->next = NULL;
        this->end = x->last;
    }
    else if (x->last == NULL)
    {
        x->next->last = NULL;
        this->start = x->next;
    }
    else
    {
        x->last->next = x->next;
        x->next->last = x->last;
    }
    delete x;
    return true;
}

// void PartitionRefinement::refine(vector<int> arr)
// {
//     vector<PartitionElement *> refinedEls;
//     for (auto x : arr)
//     {
//         PartitionElement *el = at[x];
//         if (el->refined.size() == 0)
//             refinedEls.push_back(el);
//         el->refined.insert(x);
//         el->elements.erase(x);
//     }
//     for (auto el : refinedEls)
//     {
//         PartitionElement *newel = new PartitionElement();
//         newel->elements = unordered_set<int>(el->refined.begin(), el->refined.end());
//         el->refined.clear();
//         newel->last = this->end;
//         newel->last->next = newel;
//         this->end = newel;
//         for (auto x : newel->elements)
//             at[x] = newel;

//         if (el->elements.size() == 0)
//         {
//             if (el->next == NULL)
//             {
//                 el->last->next = NULL;
//                 this->end = el->last;
//             }
//             else if (el->last == NULL)
//             {
//                 el->next->last = NULL;
//                 this->start = el->next;
//             }
//             else
//             {
//                 el->last->next = el->next;
//                 el->next->last = el->last;
//             }
//             delete el;
//         }
//     }
// }

void PartitionRefinement::print()
{
    auto el = this->start;
    while (el != NULL)
    {
        for (auto x : el->elements)
        {
            cout << x << " ";
        }
        cout << endl;
        el = el->next;
    }
}