#ifndef __QUANTUMCHIP__
#define __QUANTUMCHIP__

#include <stdlib.h>
#include <string>
#include <vector>
#include <unordered_map>

class Quantum_chip
{
    public:
        Quantum_chip(std::string, char *);
        //~Quantum_chip();
        //load coupling map from file
        int load(const char *);
        void calc_dist_m();
        std::string name; //chip name
        std::vector<std::string> gate_set;
        std::string double_gate;
        int qubit_number;
        std::vector<std::vector<int>> coupling_map; // matrix of coupling map, coupling_map[x][y] != 0 means qubit_x is able to control qubit_y and the value is gate time
        std::unordered_map<std::string,int> gate_id;
        std::vector<std::vector<int>> single_gate_time; // matrix of single qubit gate time, single_gate_time[]
        std::vector<std::vector<int>> dist;
        //std::vector<std::vector<std::vector<int>>> shortest_path;
        std::vector<std::vector<int>> steps_dist;
        std::vector<std::vector<int>> fidelity_dist;
        std::vector<std::vector<int>> neibors;
        int min_cx,max_cx,max_T1;

        std::vector<int> qubit_T1;
        std::vector<int> qubit_T2;
        std::vector<std::vector<int>> single_gate_fidelity;
        std::vector<std::vector<int>> double_gate_fidelity;
        
};
#endif