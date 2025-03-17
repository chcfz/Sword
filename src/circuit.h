#ifndef __CIRCUIT__
#define __CIRCUIT__
#include <vector>
#include "defines.h"
using std::vector;




class Gate;

class Node
{
    public:
        int id;
        int front_cnt;
        vector<Node *> successors;
        Gate *cx_gate;
        vector<Gate *> control_h_gate_list;
        vector<Gate *> target_h_gate_list;
        int control_h_to_do_index;
        int target_h_to_do_index;
        int layer;
        int has_done;
        Node(int i):id(i),front_cnt(0),cx_gate(nullptr),control_h_to_do_index(0),target_h_to_do_index(0),layer(0),has_done(false){}
        void add_successors(Node *);
};


class Gate
{
    public:
        Gate(Gate_type g_t, int i, int c, int t):
            type(g_t),id(i), control(c),target(t),owner(nullptr){}
        Gate_type type;
        int id;
        int control;
        int target;
        Node* owner;
};


class Circuit
{
    public:
        Circuit();
        Circuit(char*);
        Circuit& operator=(const Circuit &);
        ~Circuit();
        bool load(char *);
        int qubit_number;
        int type;
        vector<Node*> node_list;
        vector<Gate*> gate_list;
        vector<vector<Gate*>> left_h_gate_each_bit;
        vector<vector<int>> logical_bit_node_list;
};
#endif