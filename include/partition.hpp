#pragma once

#include "util.hpp"
#include "instance.hpp"

struct PartitionElement{
    vector<u_int32_t> elements;
    vector<u_int32_t> needrefineby;
    PartitionElement* next = NULL;
    PartitionElement* last = NULL;
};

struct PartitionRefinement{
    PartitionElement* start;
    PartitionElement* end;
    deque<PartitionElement*> nonemptyrefineby;

    PartitionRefinement(vector<u_int32_t> V, u_int32_t v);
    ~PartitionRefinement() {
        PartitionElement* el = start;
        while (el != NULL) {
            auto nextel = el->next;
            delete el;
            el = nextel;
        }
    }

    bool dorefine(OLCMInstance& instance);

    void refine(vector<int> arr);
    void print();
};