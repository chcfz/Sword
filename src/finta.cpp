#include <iostream>
#include <string.h>
#include "quantum_chip.h"
#include "circuit.h"
#include "token_swap.cpp"
#include <iterator>
#include <unordered_set>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <ctime>
#include <deque>
#define M_E 2.718281828459045
#include <cmath>
#include <regex>
#include <queue>
//const static float EXTEND_RATE = 0.5;
//const static int LOOK_AHEAD_DEPTH = 10;
float DECAY_WEIGHT = 1.0;
float DECAY_WEIGHT_FAST = 0.1;
float SINGLE_GROUP_TH = 0.25;
float SWAP_PRIOR_TH = 0.4;
float NODE_PRIOR_TH = 0.2;
int K = 10;
int NEAREND_LAYER = 10;
int ROUTETOKEN_DELTA = 1;

class TopKSmallest {
public:
    int k;
    std::priority_queue<Gate*> max_heap;  // 最大堆，堆顶是最大的元素
    TopKSmallest(int k) : k(k) {}

    // 插入数据，仅维护前K小的元素
    void add(Gate* gate) {
        if (max_heap.size() < k) {
            max_heap.push(gate);
        } else if (gate->score < max_heap.top()->score) {
            max_heap.pop();      // 移除最大的元素
            max_heap.push(gate);  // 插入更小的元素
        }
    }

    std::vector<Gate*> get() {
        std::vector<Gate*> result;
        while (!max_heap.empty()) {
            result.push_back(max_heap.top());
            max_heap.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }
};



std::vector<Circuit *> load_circuits(char *file_path) {
    std::ifstream file(file_path); // 打开文件
    if (!file.is_open()) {
        std::cout << "Failed to open file." << std::endl;
        exit(0);
    }

    std::vector<Circuit *> circuits;
    Circuit *circuit = new Circuit();
    circuit->type = Main;
    std::string line;
    std::smatch matches;
    while (std::getline(file, line)) {
        std::regex qreg_pattern(R"(\wreg \w\[(\d+)\];)");
        if (std::regex_match(line, matches, qreg_pattern)) {
            if (line[0] == 'q')
                circuit->qubit_number = std::stoi(matches[1].str());
            else {
                std::getline(file, line);
                break;
            }
        } else if (circuit->qubit_number > 0) break;
    }

    //this->qubit_gate_list = vector<vector<Gate *>>(qubit_number, vector<Gate *>());

    std::regex gate_pattern(R"(.+?\s+q\[(\d+)\](?:,\s*q\[(\d+)\])?;)");
    std::regex control_start_pattern(R"(\s*(\w+)\s*\(.+?\)\s*\{)");
    std::regex control_end_pattern(R"(\s*?})");
    Gate* gate;
    Node* node;
    Node* qbit_last_node[circuit->qubit_number] = {nullptr};
    circuit->logical_bit_node_list = vector<vector<int>>(circuit->qubit_number, vector<int>());
    vector<vector<Gate *>> tmp_h_gate_list = vector<vector<Gate *>>(circuit->qubit_number, vector<Gate *>());
    do { // 逐行读取文件
        if (std::regex_match(line, matches, gate_pattern)) { // 匹配每一行
            //std::cout << "Match found: " << matches[0] << std::endl;
            
            // 将捕获的数字字符串转换为整数
            int q1 = std::stoi(matches[1].str());

            // 检查并处理第二个捕获组
            if (matches[2].matched) {
                int q2 = std::stoi(matches[2].str());
                gate = new Gate(CX, circuit->gate_list.size(), q1, q2);
                node = new Node(circuit->node_list.size());
                node->cx_gate = gate;
                gate->owner = node;
                node->control_h_gate_list.swap(tmp_h_gate_list[q1]);
                node->target_h_gate_list.swap(tmp_h_gate_list[q2]);
                if (qbit_last_node[q1]) {
                    qbit_last_node[q1]->add_successors(node);
                }
                if (qbit_last_node[q2]) {
                    qbit_last_node[q2]->add_successors(node);
                }
                qbit_last_node[q1] = node;
                qbit_last_node[q2] = node;
                circuit->node_list.push_back(node);
                circuit->logical_bit_node_list[q1].push_back(node->id);
                circuit->logical_bit_node_list[q2].push_back(node->id);
                //std::cout << "CX " << q1 << q2 << std::endl;
                
            } else {
                //std::cout<< "H " << q1 << std::endl;
                if (line[0] == 't')
                    gate = new Gate(T, circuit->gate_list.size(), q1, -1);
                else
                    gate = new Gate(H, circuit->gate_list.size(), q1, -1);
                tmp_h_gate_list[q1].push_back(gate);
            }
            circuit->gate_list.push_back(gate);
            //this->qubit_gate_list[gate->control].push_back(gate);
            //if (gate->target != -1) this->qubit_gate_list[gate->target].push_back(gate);
        } else if (std::regex_match(line, matches, control_start_pattern)) {
            circuit->left_h_gate_each_bit.swap(tmp_h_gate_list);
            tmp_h_gate_list = vector<vector<Gate *>>(circuit->qubit_number, vector<Gate *>());
            circuits.push_back(circuit);
            circuit = new Circuit();
            circuit->qubit_number = circuits[0]->qubit_number;
            circuit->logical_bit_node_list = vector<vector<int>>(circuit->qubit_number, vector<int>());
            if (line[0] == 'i') {
                circuit->type = If;
            } else {
                circuit->type = While;
            }

            for (int i = 0; i<circuit->qubit_number; i++) {
                qbit_last_node[i] = nullptr;
            }

        } else if (std::regex_match(line, matches, control_end_pattern) ) {
            circuit->left_h_gate_each_bit.swap(tmp_h_gate_list);
            tmp_h_gate_list = vector<vector<Gate *>>(circuit->qubit_number, vector<Gate *>());
            circuits.push_back(circuit);
            circuit = new Circuit();
            circuit->qubit_number = circuits[0]->qubit_number;
            circuit->logical_bit_node_list = vector<vector<int>>(circuit->qubit_number, vector<int>());
            circuit->type = Main;

            for (int i = 0; i<circuit->qubit_number; i++) {
                qbit_last_node[i] = nullptr;
            }
        }
    } while (std::getline(file, line));
    circuit->left_h_gate_each_bit.swap(tmp_h_gate_list);
    circuits.push_back(circuit);
    for (auto c:circuits) {
        c->get_depth();
    }
    file.close(); // 关闭文件
    return circuits;
}

void dump_result(char * result_path, std::deque<Gate> &output, int qubit_number, std::vector<double> run_time_list, std::vector<int> max_decay_list) {
    // 获取文件所在的文件夹路径
    std::filesystem::path dir = std::filesystem::path(result_path).parent_path();

    // 如果文件夹不存在，则创建文件夹
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }

    // 以覆盖的方式写入文件
    std::ofstream file(result_path, std::ios::out | std::ios::trunc);

    if (file.is_open()) {
        file << "//";
        for (double run_time:run_time_list) {
            file << run_time << " ";
        }
        file << "\n";
        file << "//";
        for (int max_decay:max_decay_list) {
            file << max_decay << " ";
        }
        file << "\n";
        {
        std::ostringstream line;
        line << "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q["
             << qubit_number
             << "];\ncreg c["
             << qubit_number
             << "];\n";
        file << line.str();
        }
        
        for(auto &gate:output) {
            if (gate.id == -1) continue;
            std::ostringstream line;
            line << Gate_view[gate.type]
                 << "q[" << gate.control << "]";
            if (gate.type != H) {
                line << ",q[" << gate.target << "]";
            }
            line << ";\n";
            file << line.str();
        }

        file.close();
        //std::cout << "circuit run time: " << max_decay << std::endl;
    } else {
        std::cerr << "Failed to open the file." << std::endl;
    }
}

class Mapper {
    public:
        Quantum_chip *chip;
        Circuit *circuit;
        vector<int> decay;
        vector<int> real_decay;
        int max_decay;
        int min_decay;
        int min_decay_qubit;
        vector<Node*> front_layer;
        vector<Node*> extended_layer;
        vector<int> l2p_layout; //logical->physcial
        vector<int> p2l_layout;
        vector<int> last_swap_pair;
        vector<vector<Gate*>> last_swap_gate;
        vector<int> logical_qubit_in_front_layer_node_index; // 用来标注logical qubit对应的node在front layer中的位置
        vector<int> logical_qubit_in_extended_layer_node_index; // 用来标注logical qubit对应的node在extended layer中的位置
        std::deque<Gate> *output;
        bool routing_fail;

        Mapper(Quantum_chip *chip,
               Circuit *circuit,
               vector<int> initial_layout,
               std::deque<Gate> *output)
        {
            this->chip = chip;
            this->circuit = circuit;
            this->output = output;
            if (initial_layout.size() != 0) {
                this->l2p_layout = vector<int>(circuit->qubit_number, -1);
                for(int i = 0; i<circuit->qubit_number; i++)
                    this->l2p_layout[i] = initial_layout[i];
            }
            else {
                this->l2p_layout = vector<int>(circuit->qubit_number,-1);
                for(int i = 0; i<circuit->qubit_number; i++)
                    this->l2p_layout[i] = i;
            }

            this->p2l_layout = vector<int>(chip->qubit_number, -1);
            for (size_t i = 0; i<this->l2p_layout.size(); i++)
                this->p2l_layout[l2p_layout[i]] = i;
            for (auto node:circuit->node_list) {
                if (node->front_cnt == 0) {
                    this->front_layer.push_back(node);
                }                   
                else if (node->layer == 1 ) {
                    this->extended_layer.push_back(node);
                }
            }

            this->decay = vector<int>(chip->qubit_number ,0);
            this->real_decay = vector<int>(chip->qubit_number, 0);
            this->max_decay = 0;
            this->min_decay = 0;
            this->min_decay_qubit = 0;
            this->last_swap_pair = vector<int>(chip->qubit_number, -1);
            this->last_swap_gate = vector<vector<Gate*>>(chip->qubit_number);
            this->logical_qubit_in_front_layer_node_index = vector<int>(circuit->qubit_number, -1);
            this->logical_qubit_in_extended_layer_node_index = vector<int>(circuit->qubit_number, -1);
            this->routing_fail = false;
        }

        void reload(
            Circuit *circuit
            ) {
            this->circuit = circuit;

            for (auto node:circuit->node_list) {
                if (node->front_cnt == 0) {
                    this->front_layer.push_back(node);
                }                   
                else if (node->layer == 1 ) {
                    this->extended_layer.push_back(node);
                }
            }
            this->last_swap_pair = vector<int>(chip->qubit_number, -1);
            this->last_swap_gate = vector<vector<Gate*>>(chip->qubit_number);
            this->routing_fail = false;
        }

        vector<vector<Gate>> get_operations(vector<int> &key_qubit){
            //[bit][op]
            vector<vector<Gate>> operations;
            vector<int> has_process(this->chip->qubit_number, 0);

            auto get_swaps = [&](vector<Gate> &tmp_ops, int phy_bit) {
                for (int neib:this->chip->neibors[phy_bit]) {
                    if (has_process[neib])
                        continue;
                    tmp_ops.push_back(Gate(SWAP, 0, phy_bit, neib));
                }
                has_process[phy_bit] = 1;
                                
                if (tmp_ops.size()) {
                    operations.push_back(std::move(tmp_ops));
                    key_qubit.push_back(phy_bit);
                }

            };

            for (Node* node:this->front_layer) {
                int phy_bit_ctl = this->l2p_layout[node->cx_gate->control];
                int phy_bit_tgt = this->l2p_layout[node->cx_gate->target];
                
                vector<Gate> tmp_ops1,tmp_ops2;
                //判断是否可以执行该CX
                get_swaps(tmp_ops1, phy_bit_ctl);
                get_swaps(tmp_ops2, phy_bit_tgt);
                
            }
            return operations;
            
        }

        vector<vector<Gate>> get_operations_sort(vector<int> &key_qubit){
            //[bit][op]
            vector<vector<Gate>> operations;
            vector<int> has_process(this->chip->qubit_number, 0);

            auto get_swaps = [&](vector<Gate> &tmp_ops, int phy_bit, int pair_phy_bit) {
                for (int neib:this->chip->neibors[phy_bit]) {
                    if (has_process[neib])
                        continue;
                    tmp_ops.push_back(Gate(SWAP, 0, phy_bit, neib));
                }
                has_process[phy_bit] = 1;
                                
                if (tmp_ops.size()) {

                    std::sort(tmp_ops.begin(), tmp_ops.end(), [this, pair_phy_bit](Gate &i, Gate &j) {
                        return this->chip->dist[i.target][pair_phy_bit] < this->chip->dist[j.target][pair_phy_bit];
                    });

                    operations.push_back(std::move(tmp_ops));
                    key_qubit.push_back(phy_bit);
                }

            };
            int to_check_node_cnt = int(this->front_layer.size()*NODE_PRIOR_TH+1);

            std::sort(this->front_layer.begin(), this->front_layer.end(), [this](Node* a, Node* b) {
                return (this->logical_qubit_in_extended_layer_node_index[a->cx_gate->control] != -1 +
                       this->logical_qubit_in_extended_layer_node_index[a->cx_gate->target] != -1)
                       >
                       (this->logical_qubit_in_front_layer_node_index[b->cx_gate->control] != -1 +
                       this->logical_qubit_in_front_layer_node_index[b->cx_gate->target] != -1);
            });

            Node* node;
            for (int i = 0; i<to_check_node_cnt; i++) {
                node = this->front_layer[i];
                int phy_bit_ctl = this->l2p_layout[node->cx_gate->control];
                int phy_bit_tgt = this->l2p_layout[node->cx_gate->target];
                
                vector<Gate> tmp_ops1,tmp_ops2;
                //判断是否可以执行该CX
                get_swaps(tmp_ops1, phy_bit_ctl, phy_bit_tgt);
                get_swaps(tmp_ops2, phy_bit_tgt, phy_bit_ctl);
                
            }
            return operations;
            
        }

        vector<vector<Gate>> get_operations_sort_token(vector<int> &key_qubit, const vector<int> &target_l2p){
            //[bit][op]
            vector<vector<Gate>> operations;
            vector<int> has_process(this->chip->qubit_number, 0);

            auto get_swaps = [&](vector<Gate> &tmp_ops, int phy_bit) {
                for (int neib:this->chip->neibors[phy_bit]) {
                    if (has_process[neib])
                        continue;
                    tmp_ops.push_back(Gate(SWAP, 0, phy_bit, neib));
                }
                has_process[phy_bit] = 1;
                                
                if (tmp_ops.size()) {

                    std::sort(tmp_ops.begin(), tmp_ops.end(), [this, target_l2p](Gate &i, Gate &j) {
                        int phy1 = i.control;
                        int phy2 = i.target;
                        int logical_q1 = this->p2l_layout[phy1];
                        int logical_q2 = this->p2l_layout[phy2];
                        int lq1_target_phyq = target_l2p[logical_q1];
                        int lq2_target_phyq = target_l2p[logical_q2];
                        int delta_dist_i = this->chip->dist[lq1_target_phyq][phy2] - this->chip->dist[lq1_target_phyq][phy1]
                                    +this->chip->dist[lq2_target_phyq][phy1] - this->chip->dist[lq2_target_phyq][phy2];
                        phy1 = j.control;
                        phy2 = j.target;
                        logical_q1 = this->p2l_layout[phy1];
                        logical_q2 = this->p2l_layout[phy2];
                        lq1_target_phyq = target_l2p[logical_q1];
                        lq2_target_phyq = target_l2p[logical_q2];
                        int delta_dist_j = this->chip->dist[lq1_target_phyq][phy2] - this->chip->dist[lq1_target_phyq][phy1]
                                    +this->chip->dist[lq2_target_phyq][phy1] - this->chip->dist[lq2_target_phyq][phy2];
                        return delta_dist_i < delta_dist_j;
                    });

                    operations.push_back(std::move(tmp_ops));
                    key_qubit.push_back(phy_bit);
                }

            };
            int to_check_node_cnt = int(this->front_layer.size()*NODE_PRIOR_TH+1);
            Node* node;
            for (int i = 0; i<to_check_node_cnt; i++) {
                node = this->front_layer[i];
                int phy_bit_ctl = this->l2p_layout[node->cx_gate->control];
                int phy_bit_tgt = this->l2p_layout[node->cx_gate->target];
                
                vector<Gate> tmp_ops1,tmp_ops2;
                //判断是否可以执行该CX
                get_swaps(tmp_ops1, phy_bit_ctl);
                get_swaps(tmp_ops2, phy_bit_tgt);
                
            }
            return operations;
            
        }

        void route() {
            this->update_front_layer();
            float routing_flag_rate = 0;// the range is (0,0.5]   0.5->group choose  0->single choose
            while (this->front_layer.size() > 0) {
                vector<int> key_qbit;
                auto ready_ops = this->get_operations(key_qbit);     
                routing_flag_rate = 1.0*this->front_layer.size()/this->circuit->qubit_number;
                //vector<const Gate*> min_score_comb;
                if (routing_flag_rate > SINGLE_GROUP_TH)
                    this->routing_group(key_qbit, ready_ops);
                else
                    this->routing_single(ready_ops);
                //this->apply_comb(min_score_comb);
                this->check_front_layer_change();
            }

            for (size_t i = 0; i<circuit->left_h_gate_each_bit.size(); i++){
                for (Gate *gate:circuit->left_h_gate_each_bit[i]) {
                    output->push_back(Gate(H, gate->id, l2p_layout[gate->control],-1));
                }
            }
        }

        void routing_group(vector<int> &key_qbit, vector<vector<Gate>> &ready_ops) {
            static vector<int> bit_has_added(this->chip->qubit_number, 0);
            bit_has_added.assign(this->chip->qubit_number, 0);
            Gate *min_score_gate;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist;
            int now_extend_dist;
            while (1)
            {
                min_score_gate = nullptr;
                now_front_dist = this->front_dist(this->chip->dist);
                now_extend_dist = this->extend_dist(this->chip->dist);
                for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                    if (bit_has_added[key_qbit[i]]) continue;
                    for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()); j++){
                        if (bit_has_added[ready_ops[i][j].control] || bit_has_added[ready_ops[i][j].target]) continue;
                        //min_score_comb.push_back(&ready_ops[i][j]);
                        tmp_score = this->score_single(&ready_ops[i][j], now_front_dist, now_extend_dist);
                        if (min_score>tmp_score) {
                            min_score = tmp_score;
                            //min_reduce_delta = tmp_delta;
                            min_score_gate = &ready_ops[i][j];
                        }
                        //min_score_comb.pop_back();
                    }
                }
                if (!min_score_gate) break;
                //min_score_comb.push_back(min_score_gate);
                apply_swap(min_score_gate);
                bit_has_added[min_score_gate->control] = 1;
                bit_has_added[min_score_gate->target] = 1;
                //break;
            }
        }

        void routing_single(vector<vector<Gate>> &ready_ops) {
            Gate *min_score_gate = nullptr;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist = this->front_dist(this->chip->dist);
            int now_extend_dist = this->extend_dist(this->chip->dist);
            for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()); j++){
                    tmp_score = this->score_single(&ready_ops[i][j], now_front_dist, now_extend_dist);
                    if (min_score>tmp_score) {
                        min_score = tmp_score;
                        min_score_gate = &ready_ops[i][j];
                    }
                }
            }
            if (min_score_gate)
               apply_swap(min_score_gate);
        }

        float score_single(Gate* op, int now_front_dist, int now_extend_dist) {
            int phy_q1,phy_q2;
            int part_max_decay = 0;
            phy_q1 = op->control;
            phy_q2 = op->target;
            part_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
            
            int front_score = now_front_dist + this->dist_change_after_swap(op, logical_qubit_in_front_layer_node_index, front_layer);
            int extend_score = now_extend_dist + this->dist_change_after_swap(op, logical_qubit_in_extended_layer_node_index, extended_layer);
            float rate = 1.0;
            if (this->max_decay != 0) {
                if (min_decay < max_decay)
                    rate += DECAY_WEIGHT*(part_max_decay-min_decay)/(max_decay-min_decay);
                else
                    rate += DECAY_WEIGHT*part_max_decay/max_decay;
            }
            if (extended_layer.size() == 0) {
                return  rate * (1.0*front_score/front_layer.size());
            } else {
                return  rate * (1.0*front_score/front_layer.size() + 0.1*extend_score/extended_layer.size());
            }
        }

        void route_fast(const vector<int> &target_l2p) {
            this->update_front_layer();
            float routing_flag_rate = 0;// the range is (0,0.5]   0.5->group choose  0->single choose
            while (this->front_layer.size() > 0) {
                vector<int> key_qbit;
                auto ready_ops = this->get_operations_sort(key_qbit);
                routing_flag_rate = 1.0*this->front_layer.size()/this->circuit->qubit_number;
                //vector<const Gate*> min_score_comb;
                if (routing_flag_rate > SINGLE_GROUP_TH)
                    this->routing_group_fast(key_qbit, ready_ops);
                else
                    this->routing_single_fast(ready_ops);
                //this->apply_comb(min_score_comb);
                this->check_front_layer_change();
                if (this->check_nearing_end()) {
                    break;
                }
            }

            while (this->front_layer.size() > 0) {
                vector<int> key_qbit;
                auto ready_ops = this->get_operations_sort(key_qbit);
                routing_flag_rate = 1.0*this->front_layer.size()/this->circuit->qubit_number;
                //vector<const Gate*> min_score_comb;
                if (routing_flag_rate > SINGLE_GROUP_TH)
                    this->routing_group_fast_token(key_qbit, ready_ops, target_l2p);
                else
                    this->routing_single_fast_token(ready_ops, target_l2p);
                //this->apply_comb(min_score_comb);
                this->check_front_layer_change();
            }

            for (size_t i = 0; i<circuit->left_h_gate_each_bit.size(); i++){
                for (Gate *gate:circuit->left_h_gate_each_bit[i]) {
                    output->push_back(Gate(H, gate->id, l2p_layout[gate->control],-1));
                }
            }
        }
        

        void routing_group_fast(vector<int> &key_qbit, vector<vector<Gate>> &ready_ops) {
            static vector<int> bit_has_added(this->chip->qubit_number, 0);
            bit_has_added.assign(this->chip->qubit_number, 0);
            Gate *min_score_gate;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist;
            int now_extend_dist;
            while (1)
            {
                min_score_gate = nullptr;
                now_front_dist = this->front_dist(this->chip->dist);
                now_extend_dist = this->extend_dist(this->chip->dist);
                for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                    if (bit_has_added[key_qbit[i]]) continue;
                    for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                        if (bit_has_added[ready_ops[i][j].control] || bit_has_added[ready_ops[i][j].target]) continue;
                        //min_score_comb.push_back(&ready_ops[i][j]);
                        tmp_score = this->score_single_fast(&ready_ops[i][j], now_front_dist, now_extend_dist);
                        if (min_score>tmp_score) {
                            min_score = tmp_score;
                            //min_reduce_delta = tmp_delta;
                            min_score_gate = &ready_ops[i][j];
                        }
                        //min_score_comb.pop_back();
                    }
                }
                if (!min_score_gate) break;
                //min_score_comb.push_back(min_score_gate);
                apply_swap(min_score_gate);
                bit_has_added[min_score_gate->control] = 1;
                bit_has_added[min_score_gate->target] = 1;
                //break;
            }
        }

        void routing_single_fast(vector<vector<Gate>> &ready_ops) {
            Gate *min_score_gate = nullptr;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist = this->front_dist(this->chip->dist);
            int now_extend_dist = this->extend_dist(this->chip->dist);
            for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                    tmp_score = this->score_single_fast(&ready_ops[i][j], now_front_dist, now_extend_dist);
                    if (min_score>tmp_score) {
                        min_score = tmp_score;
                        min_score_gate = &ready_ops[i][j];
                    }
                }
            }
            if (min_score_gate)
               apply_swap(min_score_gate);
        }

        float score_single_fast(Gate* op, int now_front_dist, int now_extend_dist) {
            int phy_q1,phy_q2;
            int part_max_decay = 0;
            phy_q1 = op->control;
            phy_q2 = op->target;
            part_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
            
            int front_score = now_front_dist + this->dist_change_after_swap(op, logical_qubit_in_front_layer_node_index, front_layer);
            int extend_score = now_extend_dist + this->dist_change_after_swap(op, logical_qubit_in_extended_layer_node_index, extended_layer);
            float rate = 1.0;
            if (this->max_decay != 0) {
                if (min_decay < max_decay)
                    rate += DECAY_WEIGHT_FAST*(part_max_decay-min_decay)/(max_decay-min_decay);
                else
                    rate += DECAY_WEIGHT_FAST*part_max_decay/max_decay;
            }
            if (extended_layer.size() == 0) {
                return  rate * (1.0*front_score/front_layer.size());
            } else {
                return  rate * (1.0*front_score/front_layer.size() + 0.1*extend_score/extended_layer.size());
            }
        }


        
        void routing_group_fast_token(vector<int> &key_qbit, vector<vector<Gate>> &ready_ops, const vector<int> &target_l2p) {
            static vector<int> bit_has_added(this->chip->qubit_number, 0);
            bit_has_added.assign(this->chip->qubit_number, 0);
            Gate *min_score_gate;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist;
            int now_extend_dist;
            while (1)
            {
                min_score_gate = nullptr;
                now_front_dist = this->front_dist(this->chip->dist);
                now_extend_dist = this->extend_dist(this->chip->dist);

                TopKSmallest top_k_smallest(K);

                for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                    if (bit_has_added[key_qbit[i]]) continue;
                    for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                        if (bit_has_added[ready_ops[i][j].control] || bit_has_added[ready_ops[i][j].target]) continue;
                        //min_score_comb.push_back(&ready_ops[i][j]);
                        tmp_score = this->score_single_fast(&ready_ops[i][j], now_front_dist, now_extend_dist);
                        ready_ops[i][j].score = tmp_score;
                        top_k_smallest.add(&ready_ops[i][j]);
                        if (min_score>tmp_score) {
                            min_score = tmp_score;
                            min_score_gate = &ready_ops[i][j];
                        }
                    }
                }
                if (!min_score_gate) break;
                
                int phy1, phy2;
                int logical_q1, logical_q2;
                int lq1_target_phyq, lq2_target_phyq;
                int delta_dist = 0;
                int min_delta_dist = INT32_MAX;
                for (Gate* gate:top_k_smallest.get()) {
                    if (gate->score - min_score > ROUTETOKEN_DELTA) break;
                    phy1 = gate->control;
                    phy2 = gate->target;
                    logical_q1 = this->p2l_layout[phy1];
                    logical_q2 = this->p2l_layout[phy2];
                    lq1_target_phyq = target_l2p[logical_q1];
                    lq2_target_phyq = target_l2p[logical_q2];
                    delta_dist = chip->dist[lq1_target_phyq][phy2] - chip->dist[lq1_target_phyq][phy1]
                                +chip->dist[lq2_target_phyq][phy1] - chip->dist[lq2_target_phyq][phy2];
                    if (min_delta_dist > delta_dist) {
                        min_delta_dist = delta_dist;
                        min_score_gate = gate;
                    }
                }

                apply_swap(min_score_gate);
                bit_has_added[min_score_gate->control] = 1;
                bit_has_added[min_score_gate->target] = 1;
            }
        }

        void routing_single_fast_token(vector<vector<Gate>> &ready_ops, const vector<int> &target_l2p) {
            Gate *min_score_gate = nullptr;
            float min_score = FLT_MAX;
            float tmp_score;
            int now_front_dist = this->front_dist(this->chip->dist);
            int now_extend_dist = this->extend_dist(this->chip->dist);
            TopKSmallest top_k_smallest(K);
            for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                    tmp_score = this->score_single_fast(&ready_ops[i][j], now_front_dist, now_extend_dist);
                    ready_ops[i][j].score = tmp_score;
                    top_k_smallest.add(&ready_ops[i][j]);
                    if (min_score>tmp_score) {
                        min_score = tmp_score;
                        min_score_gate = &ready_ops[i][j];
                    }
                }
            }
            if (min_score_gate) {
                int phy1, phy2;
                int logical_q1, logical_q2;
                int lq1_target_phyq, lq2_target_phyq;
                int delta_dist = 0;
                int min_delta_dist = INT32_MAX;
                for (Gate* gate:top_k_smallest.get()) {
                    if (gate->score - min_score > ROUTETOKEN_DELTA) break;
                    phy1 = gate->control;
                    phy2 = gate->target;
                    logical_q1 = this->p2l_layout[phy1];
                    logical_q2 = this->p2l_layout[phy2];
                    lq1_target_phyq = target_l2p[logical_q1];
                    lq2_target_phyq = target_l2p[logical_q2];
                    delta_dist = chip->dist[lq1_target_phyq][phy2] - chip->dist[lq1_target_phyq][phy1]
                                +chip->dist[lq2_target_phyq][phy1] - chip->dist[lq2_target_phyq][phy2];
                    if (min_delta_dist > delta_dist) {
                        min_delta_dist = delta_dist;
                        min_score_gate = gate;
                    }
                }
                apply_swap(min_score_gate);
            }
        }


        int dist_change_after_swap(Gate* op, vector<int> &logical_qubit_in_layer_node_index, vector<Node*> &layer) { 
            int logical_q1 = this->p2l_layout[op->control];
            int logical_q2 = this->p2l_layout[op->target];
            int node_id1 = logical_qubit_in_layer_node_index[logical_q1];
            int node_id2 = logical_qubit_in_layer_node_index[logical_q2];
            Node *node;

            int target_phy_qubit;
            int delta_dist = 0;
            if (node_id1 != -1) {
                node = layer[node_id1];
                if (node->cx_gate->control == logical_q1) { 
                    target_phy_qubit = l2p_layout[node->cx_gate->target];
                } else {
                    target_phy_qubit = l2p_layout[node->cx_gate->control];
                }
                delta_dist += chip->dist[target_phy_qubit][op->target] - chip->dist[target_phy_qubit][op->control];
            }
            if (node_id2 != -1 && node_id2 != node_id1) {
                node = layer[node_id2];
                if (node->cx_gate->control == logical_q2) { 
                    target_phy_qubit = l2p_layout[node->cx_gate->target];
                } else {
                    target_phy_qubit = l2p_layout[node->cx_gate->control];
                }
                delta_dist += chip->dist[target_phy_qubit][op->control] - chip->dist[target_phy_qubit][op->target];
            }
            return delta_dist;
        }
        
        int front_dist(vector<vector<int>> &dist_m) {
            int phy_q1,phy_q2;
            int total_end_time = 0;
            Gate* gate;
            for (auto node:this->front_layer) {
                gate = node->cx_gate;
                phy_q1 = this->l2p_layout[gate->control];
                phy_q2 = this->l2p_layout[gate->target];

                total_end_time += dist_m[phy_q1][phy_q2];

            }
            return total_end_time;
        }

        int extend_dist(vector<vector<int>> &dist_m) {
            if (extended_layer.size()==0) return 0;
            int sum_dist = 0;
            int q1,q2;
            for (auto node:this->extended_layer) {
                q1 = this->l2p_layout[node->cx_gate->control];
                q2 = this->l2p_layout[node->cx_gate->target];
                sum_dist += dist_m[q1][q2];
            }
            return sum_dist;
        }

        void swap(const Gate *gate) {
            if (gate->type != SWAP) return;
            int phy_q1 = gate->control;
            int phy_q2 = gate->target;
            int logical_q1 = this->p2l_layout[phy_q1];
            int logical_q2 = this->p2l_layout[phy_q2];

            this->p2l_layout[phy_q1] = logical_q2;
            this->p2l_layout[phy_q2] = logical_q1;

            if(logical_q1 != -1) this->l2p_layout[logical_q1] = phy_q2;
            if(logical_q2 != -1) this->l2p_layout[logical_q2] = phy_q1;
        }

        void apply_swap(Gate* gate) {
            int phy_q1,phy_q2;
            int two_bit_max_decay;
            
            phy_q1 = gate->control;
            phy_q2 = gate->target;
            #ifdef DEBUG
            std::cout<<"swap "<<phy_q1<<" "<<phy_q2<<std::endl;
            #endif
            two_bit_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
            this->decay[phy_q1] = two_bit_max_decay + 3*this->chip->coupling_map[phy_q1][phy_q2];
            this->decay[phy_q2] = this->decay[phy_q1];
            this->max_decay = std::max(this->decay[phy_q1], this->max_decay);
            // if (phy_q1 == min_decay_qubit || phy_q2 == min_decay_qubit) {
            //     for (size_t i = 0; i<this->decay.size(); i++) {
            //         if (decay[i] < decay[min_decay_qubit]) min_decay_qubit = i;
            //     }
            //     min_decay = decay[min_decay_qubit];
            // }
            
            swap(gate);
            if (last_swap_pair[phy_q1] == phy_q2 && last_swap_pair[phy_q2] == phy_q1) {
                last_swap_gate[phy_q1].back()->id = -1;
                last_swap_gate[phy_q1].pop_back();
                last_swap_gate[phy_q2].pop_back();

                real_decay[phy_q1] -= 3*chip->coupling_map[phy_q1][phy_q2];
                real_decay[phy_q2] = real_decay[phy_q1];

                Gate* swap_gate;
                if (! last_swap_gate[phy_q1].empty()) {
                    swap_gate = last_swap_gate[phy_q1].back();
                    last_swap_pair[phy_q1] = phy_q1 != swap_gate->control?swap_gate->control:swap_gate->target;
                } else {
                    last_swap_pair[phy_q1] = -1;
                }
                if (! last_swap_gate[phy_q2].empty()) {
                    swap_gate = last_swap_gate[phy_q2].back();
                    last_swap_pair[phy_q2] = phy_q2 != swap_gate->control?swap_gate->control:swap_gate->target;
                } else {
                    last_swap_pair[phy_q2] = -1;
                }
            } else {
                two_bit_max_decay = std::max(this->real_decay[phy_q1], this->real_decay[phy_q2]);
                this->real_decay[phy_q1] = two_bit_max_decay + 3*this->chip->coupling_map[phy_q1][phy_q2];
                this->real_decay[phy_q2] = this->real_decay[phy_q1];
                
                run_H_gate_before_op(gate->control, two_bit_max_decay);
                run_H_gate_before_op(gate->target, two_bit_max_decay);

                auto swap_gate = Gate(SWAP, 0, phy_q1, phy_q2);
                this->last_swap_pair[phy_q1] = phy_q2;
                this->last_swap_pair[phy_q2] = phy_q1;
                output->push_back(std::move(swap_gate));
                this->last_swap_gate[phy_q1].push_back(&output->back());
                this->last_swap_gate[phy_q2].push_back(&output->back());
            }



        }

        bool check_front_layer_change() {
            int phy_q1,phy_q2;
            bool front_layer_change = false;
            for (auto node:front_layer) {
                phy_q1 = l2p_layout[node->cx_gate->control];
                phy_q2 = l2p_layout[node->cx_gate->target];
                if (chip->coupling_map[phy_q1][phy_q2]) {
                    front_layer_change = true;
                    break;
                }
            }
            if (front_layer_change) this->update_front_layer();
            return front_layer_change;
        }

        int run_H_gate_before_op(int phy_qubit, int latest_end_time) {
            int logical_qubit = this->p2l_layout[phy_qubit];
            if (logical_qubit == -1) return false;
            if (this->circuit->logical_bit_node_list[logical_qubit].size() == 0) return false;
            int node_id = this->circuit->logical_bit_node_list[logical_qubit][0];
            Node* node = this->circuit->node_list[node_id];

            int* start_index;
            vector<Gate*> *start_h_list;
            if (logical_qubit == node->cx_gate->control){
                start_index = &(node->control_h_to_do_index);
                start_h_list = &(node->control_h_gate_list);
            } else {
                start_index = &(node->target_h_to_do_index);
                start_h_list = &(node->target_h_gate_list);                
            }

            int start_time = this->real_decay[phy_qubit];
            int gate_run_time;
            Gate *tmp;
            for(;*start_index<static_cast<int>(start_h_list->size());(*start_index)++) {
                tmp = (*start_h_list)[*start_index];
                gate_run_time = this->chip->single_gate_time[phy_qubit][tmp->type];
                if (start_time + gate_run_time <= latest_end_time) {
                    start_time += gate_run_time;
                    //std::cout<<"h "<<phy_qubit<< "\t|| h "<<logical_qubit<<std::endl;
                    output->push_back(Gate(tmp->type, tmp->id, phy_qubit, -1));
                } else {
                    break;
                }
            }
            //如果 start_time改变了，说明有成功执行H
            this->real_decay[phy_qubit] = start_time;
            return start_time;
        }

        void run_node(Node* node) {
            int logical_q1 = node->cx_gate->control;
            int logical_q2 = node->cx_gate->target;
            int phy_q1 = l2p_layout[logical_q1];
            int phy_q2 = l2p_layout[logical_q2];
            //std::cout<<logical_q1<<" "<<logical_q2<<" || ";
            //std::cout<<phy_q1<<" "<<phy_q2<<std::endl;
            int phy_q1_start_time = run_H_gate_before_op(phy_q1, INT32_MAX);
            int phy_q2_start_time = run_H_gate_before_op(phy_q2, INT32_MAX);
            int end_time = std::max(phy_q1_start_time, phy_q2_start_time)
                         + this->chip->coupling_map[phy_q1][phy_q2];
            
            this->last_swap_pair[phy_q1] = -1;
            this->last_swap_pair[phy_q2] = -1;
            this->last_swap_gate[phy_q1].clear();
            this->last_swap_gate[phy_q2].clear();
            this->real_decay[phy_q1] = end_time;
            this->real_decay[phy_q2] = end_time;
            auto &tmp_list1 =  this->circuit->logical_bit_node_list[logical_q1];
            auto &tmp_list2 =  this->circuit->logical_bit_node_list[logical_q2];
            tmp_list1.erase(tmp_list1.begin());
            tmp_list2.erase(tmp_list2.begin());
            
            node->has_done = true;
            //std::cout<<"CX "<<phy_q1<<" "<<phy_q2<<std::endl;
            output->push_back(Gate(CX,node->cx_gate->id, phy_q1, phy_q2));
        }

        void update_front_layer() {

            this->logical_qubit_in_front_layer_node_index.assign(this->circuit->qubit_number, -1);
            this->logical_qubit_in_extended_layer_node_index.assign(this->circuit->qubit_number, -1);
            vector<Node*> tmp;
            int phy_q1,phy_q2;
            for (size_t i =0 ; i<front_layer.size(); i++) {
                phy_q1 = l2p_layout[front_layer[i]->cx_gate->control];
                phy_q2 = l2p_layout[front_layer[i]->cx_gate->target];
                if (chip->coupling_map[phy_q1][phy_q2]) {
                    run_node(front_layer[i]);
                    for (auto node:front_layer[i]->successors) {
                        node->front_cnt--;
                        if (node->front_cnt == 0) {
                            front_layer.push_back(node);
                        }
                    }
                } else {
                    logical_qubit_in_front_layer_node_index[front_layer[i]->cx_gate->control] = tmp.size();
                    logical_qubit_in_front_layer_node_index[front_layer[i]->cx_gate->target] = tmp.size();
                    tmp.push_back(front_layer[i]);
                }
            }
            front_layer.swap(tmp);
            extended_layer.clear();
            for (auto node:front_layer)
            for (auto succ:node->successors) {
                //借用has done, 防止重复添加
                if (succ->has_done) continue;
                if (succ->front_cnt == 1 ||succ->layer == node->layer+1){
                    logical_qubit_in_extended_layer_node_index[succ->cx_gate->control] = extended_layer.size();
                    logical_qubit_in_extended_layer_node_index[succ->cx_gate->target] = extended_layer.size();

                    extended_layer.push_back(succ);
                    succ->has_done = true;
                }
            }

            for (auto node:extended_layer)
                node->has_done = false;        
        
            //同步decay 与real_decay;
            this->max_decay = 0;
            this->min_decay_qubit = 0;
            for(int i = 0; i<this->chip->qubit_number; i++) {
                decay[i] = real_decay[i];
                max_decay = std::max(max_decay, decay[i]);
                if (decay[min_decay_qubit] > decay[i]) min_decay_qubit = i;
            }
            min_decay = decay[min_decay_qubit];
        }

        bool check_nearing_end() {
            for (auto node:front_layer) {
                if (circuit->depth - node->layer <= NEAREND_LAYER) return true;
            }
            return false;
        }

};


int main(int argc, char **argv) {
    if (argc < 7) {
        return 1;
    }
    char *chip_path = argv[1];
    char *circuit_path = argv[2];
    char *result_path = argv[3];
    char *route_mode = argv[4]; // online(fast)  offline(normal)    -->  'f',  'n'
    K = atoi(argv[5]);
    NEAREND_LAYER = atoi(argv[6]);
    std::deque<Gate> output;
    std::vector<double> run_time_list;
    std::vector<int> max_decay_list;
    double run_time;
    auto circuits = load_circuits(circuit_path);
    Circuit *c = circuits[0];
    Quantum_chip *chip = new Quantum_chip("aaa", chip_path);
    vector<int> initial_layout(c->qubit_number);
    for (int i =0; i<c->qubit_number; i++)
        initial_layout[i] = i;
    Mapper mp(chip, c, initial_layout, &output);

    clock_t start_time = clock();
    mp.route();
    clock_t end_time = clock();
    run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
    run_time_list.push_back(run_time);
    max_decay_list.push_back(mp.max_decay);
    double token_swap_time = 0;
    vector last_final_mapping = mp.p2l_layout;
    vector<int> target_l2p = mp.l2p_layout;
    for (size_t i = 1; i<circuits.size(); i++) {
        mp.reload(circuits[i]);
        if (i % 2 == 0) {
            clock_t start_time = clock();
            mp.route();
            clock_t end_time = clock();
            run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
        }
        else {
            //subcircuit
            //token swap
            clock_t start_time = clock();
            if (route_mode[0] == 'f') {
                mp.route_fast(target_l2p);
            }
            else {
                mp.route();
            }
            clock_t end_time = clock();
            run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
            auto swap_sequence = TS4_new(chip->coupling_map, mp.p2l_layout, last_final_mapping, &token_swap_time);
            for (auto swap:swap_sequence) {
                auto gate = Gate(SWAP, 0, swap[0], swap[1]);
                mp.apply_swap(&gate);
            }
            run_time += token_swap_time;
            //check ts is correct
            // for (size_t i = 0; i<mp.p2l_layout.size(); i++) {
            //     if (mp.p2l_layout[i] != last_final_mapping[i]) {
            //         std::cout<<"ts is not correct"<<std::endl;
            //         exit(0);
            //     }
            // }
        }
        last_final_mapping = mp.p2l_layout;
        target_l2p = mp.l2p_layout;
        run_time_list.push_back(run_time);
        max_decay_list.push_back(mp.max_decay);
    }
    std::cout<<std::endl;
    std::cout<<result_path<<std::endl;
    dump_result(result_path, output, chip->qubit_number, run_time_list, max_decay_list);
    for (Circuit *circuit:circuits) {
        delete circuit;
    }
    delete chip;
    return 0;
}
