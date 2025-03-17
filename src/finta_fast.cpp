#include <iostream>
#include <string.h>
#include "quantum_chip.h"
#include "circuit.h"
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
#include <random>

//const static float EXTEND_RATE = 0.5;
//const static int LOOK_AHEAD_DEPTH = 10;
// float COSUME_WEIGHT = 0.1;
float DECAY_WEIGHT = 1.0;
float SINGLE_GROUP_TH = 0.25;
float SWAP_PRIOR_TH = 0.4;
float NODE_PRIOR_TH = 0.5;
bool opt_fidelity = false;


class Mapper {
    public:
        Quantum_chip *chip;
        Circuit *circuit;
        vector<int> decay;
        vector<int> real_decay;
        int max_decay;
        int free_time_delta; // = min(gate time)/2
        vector<Node*> front_layer;
        vector<Node*> extended_layer;
        vector<int> l2p_layout; //logical->physcial
        vector<int> p2l_layout;
        //vector<Node*> logical_bit_now_node;
        vector<int> last_swap_pair;
        vector<vector<Gate*>> last_swap_gate;
        std::deque<Gate> output;
        bool routing_fail;

        Mapper(Quantum_chip *chip,
               Circuit *circuit,
               vector<int> initial_layout)
        {
            this->chip = chip;
            this->circuit = circuit;
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
            this->last_swap_pair = vector<int>(chip->qubit_number, -1);
            this->last_swap_gate = vector<vector<Gate*>>(chip->qubit_number);
            this->routing_fail = false;
        }

        void dump_result(char * result_path) {
            if (routing_fail) return;
            // 获取文件所在的文件夹路径
            std::filesystem::path dir = std::filesystem::path(result_path).parent_path();

            // 如果文件夹不存在，则创建文件夹
            if (!std::filesystem::exists(dir)) {
                std::filesystem::create_directories(dir);
            }

            // 以覆盖的方式写入文件
            std::ofstream file(result_path, std::ios::out | std::ios::trunc);

            for (size_t i = 0; i<circuit->left_h_gate_each_bit.size(); i++){
                for (Gate *gate:circuit->left_h_gate_each_bit[i]) {
                    output.push_back(Gate(H, gate->id, l2p_layout[gate->control],-1));
                }
            }

            if (file.is_open()) {
                {
                std::ostringstream line;
                line << "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q["
                     << this->chip->qubit_number
                     << "];\ncreg c["
                     << this->chip->qubit_number
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

        vector<vector<Gate>> get_operations(vector<int> &key_qubit){
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
                        if (opt_fidelity) {
                            return this->chip->fidelity_dist[i.target][pair_phy_bit] < this->chip->fidelity_dist[j.target][pair_phy_bit];
                        }
                        return this->chip->dist[i.target][pair_phy_bit] < this->chip->dist[j.target][pair_phy_bit];
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
                get_swaps(tmp_ops1, phy_bit_ctl, phy_bit_tgt);
                get_swaps(tmp_ops2, phy_bit_tgt, phy_bit_ctl);
                
            }
            return operations;
            
        }

        void route() {
            this->update_front_layer();
            float routing_flag_rate = 0;// the range is (0,0.5]   0.5->group choose  0->single choose
            // std::random_device rd;
            // std::mt19937 g(rd());
            while (this->front_layer.size() > 0) {
                vector<int> key_qbit;
                auto ready_ops = this->get_operations(key_qbit);
                routing_flag_rate = 1.0*this->front_layer.size()/this->circuit->qubit_number;
                vector<const Gate*> min_score_comb;
                if (routing_flag_rate>0.5) std::cout<<routing_flag_rate<<std::endl;
                if (routing_flag_rate > SINGLE_GROUP_TH)
                    this->routing_group(key_qbit, ready_ops, min_score_comb);
                else
                    this->routing_single(ready_ops, min_score_comb);
                if (min_score_comb.size()==0) {
                    std::cout<<"error routing"<<std::endl;
                    this->routing_fail = true;
                    break;
                }
                this->apply_comb(min_score_comb);
            }
        }


        void routing_group(vector<int> &key_qbit, vector<vector<Gate>> &ready_ops, vector<const Gate*> &min_score_comb) {
            static vector<int> bit_has_added(this->chip->qubit_number, 0);
            bit_has_added.assign(this->chip->qubit_number, 0);
            Gate *min_score_gate;
            float min_score = FLT_MAX;
            float tmp_score;
            while (1)
            {
                min_score_gate = nullptr;
                //now_front_score = this->front_dist(&empty_consume_time);
                for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                    if (bit_has_added[key_qbit[i]]) continue;
                    for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                        if (bit_has_added[ready_ops[i][j].control] || bit_has_added[ready_ops[i][j].target]) continue;
                        min_score_comb.push_back(&ready_ops[i][j]);
                        if (opt_fidelity)
                            tmp_score = this->score_group_fidelity(min_score_comb);
                        else
                            tmp_score = this->score_group(min_score_comb);
                        if (min_score>tmp_score) {
                            min_score = tmp_score;
                            //min_reduce_delta = tmp_delta;
                            min_score_gate = &ready_ops[i][j];
                        }
                        min_score_comb.pop_back();
                    }
                }
                if (!min_score_gate) break;
                min_score_comb.push_back(min_score_gate);
                bit_has_added[min_score_gate->control] = 1;
                bit_has_added[min_score_gate->target] = 1;
                //break;
            }
        }

        void routing_single(vector<vector<Gate>> &ready_ops,vector<const Gate*> &min_score_comb) {
            Gate *min_score_gate = nullptr;
            float min_score = FLT_MAX;
            float tmp_score;
            for(size_t i = 0; i<static_cast<size_t>(ready_ops.size()); i++) {
                for(size_t j = 0; j<static_cast<size_t>(ready_ops[i].size()*SWAP_PRIOR_TH+1); j++){
                    if (opt_fidelity)
                        tmp_score = this->score_fidelity(&ready_ops[i][j]);
                    else
                        tmp_score = this->score_single(&ready_ops[i][j]);
                    if (min_score>tmp_score) {
                        min_score = tmp_score;
                        min_score_gate = &ready_ops[i][j];
                    }
                }
            }
            if (min_score_gate)
                min_score_comb.push_back(min_score_gate);
        }

        float score_fidelity(Gate* op) {
            static const float M_1_MAX_T1 = 1.0/(chip->max_T1);
            int phy_q1,phy_q2;
            int part_max_decay = 0;
            phy_q1 = op->control;
            phy_q2 = op->target;
            part_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
            swap(op);
            //total_consume_time += 3 * this->chip->coupling_map[phy_q1][phy_q2];
            int front_score = this->front_dist(this->chip->fidelity_dist);
            int extend_score = this->extend_dist(this->chip->fidelity_dist);
            //恢复
            swap(op);
            float rate = 1.0;
            if (this->max_decay != 0) {
                rate += DECAY_WEIGHT*part_max_decay/this->max_decay;
            }
            return  rate * (front_score + 0.1*extend_score);
        }

        float score_single(Gate* op) {
            //static const float M_1_MAX_CX_DELTA = 1.0/(chip->max_cx);
            int phy_q1,phy_q2;
            int part_max_decay = 0;
            int total_consume_time = 0;
            phy_q1 = op->control;
            phy_q2 = op->target;
            part_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
            swap(op);
            total_consume_time += 3 * this->chip->coupling_map[phy_q1][phy_q2];
            int front_score = this->front_dist(this->chip->steps_dist);
            int extend_score = this->extend_dist(this->chip->steps_dist);
            //恢复
            swap(op);
            float rate = 1.0;
            if (this->max_decay != 0) {
                rate += DECAY_WEIGHT*part_max_decay/this->max_decay;
            }
            return  rate * (front_score + 0.1*extend_score);
        }

        float score_group_fidelity(vector<const Gate*>& currentCombination) {
            static const float M_1_MAX_T1_DELTA = 1.0/(chip->max_T1);
            int phy_q1,phy_q2;
            int tmp_max_end_decay = 0;
            int part_max_decay = 0;
            //int total_consume_fidelity = 0;
            for (auto op:currentCombination) {
                phy_q1 = op->control;
                phy_q2 = op->target;
                tmp_max_end_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
                swap(op);
                tmp_max_end_decay += 3 * this->chip->coupling_map[phy_q1][phy_q2];
                part_max_decay = std::max(part_max_decay, tmp_max_end_decay);
                //total_consume_fidelity += 3 * this->chip->double_gate_fidelity[phy_q1][phy_q2];
            }

            int front_score = this->front_dist(this->chip->fidelity_dist);
            int extend_score = this->extend_dist(this->chip->fidelity_dist);
            //恢复
            for (auto op:currentCombination) {
                swap(op);
            }
            float rate = 1.0;
            if (this->max_decay != 0) {
                phy_q1 = currentCombination.back()->control;
                phy_q2 = currentCombination.back()->target;
                int max_qubit_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
                rate += DECAY_WEIGHT*max_qubit_decay/this->max_decay;
            }
            return  rate * (front_score + 0.1*extend_score);
        }


        float score_group(vector<const Gate*>& currentCombination) {
            static const float M_1_MAX_CX_DELTA = 1.0/(chip->max_cx);
            int phy_q1,phy_q2;
            int tmp_max_end_decay = 0;
            int part_max_decay = 0;
            int total_consume_time = 0;
            for (auto op:currentCombination) {
                phy_q1 = op->control;
                phy_q2 = op->target;
                tmp_max_end_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
                swap(op);
                tmp_max_end_decay += 3 * this->chip->coupling_map[phy_q1][phy_q2];
                part_max_decay = std::max(part_max_decay, tmp_max_end_decay);
                tmp_max_end_decay = 0;

                total_consume_time += 3 * this->chip->coupling_map[phy_q1][phy_q2];
            }

            int front_score = this->front_dist(this->chip->dist);
            int extend_score = this->extend_dist(this->chip->dist);
            //恢复
            for (auto op:currentCombination) {
                swap(op);
            }
            float rate = 1.0;
            if (this->max_decay != 0) {
                phy_q1 = currentCombination.back()->control;
                phy_q2 = currentCombination.back()->target;
                int max_qubit_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
                rate += DECAY_WEIGHT*max_qubit_decay/this->max_decay;
            }
            return  rate * (front_score + 0.1*extend_score);
        }
        int front_dist(vector<vector<int>> &dist_m) {
            int phy_q1,phy_q2;
            int total_end_time = 0;
            Gate* gate;
            for (auto node:this->front_layer) {
                gate = node->cx_gate;
                phy_q1 = this->l2p_layout[gate->control];
                phy_q2 = this->l2p_layout[gate->target];

                total_end_time += 1.5*dist_m[phy_q1][phy_q2];                   

            }
            return total_end_time/front_layer.size();
        }

        int extend_dist(vector<vector<int>> &dist_m) {
            if (extended_layer.size()==0) return 0;
            int sum_dist = 0;
            int q1,q2;
            for (auto node:this->extended_layer) {
                q1 = this->l2p_layout[node->cx_gate->control];
                q2 = this->l2p_layout[node->cx_gate->target];
                if (this->chip->coupling_map[q1][q2])
                    sum_dist += dist_m[q1][q2];
                else sum_dist += 1.5*dist_m[q1][q2];
            }
            return sum_dist/(extended_layer.size());
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

        bool apply_comb(vector<const Gate*>& comb) {
            static int cnt = 0;
            cnt ++;
            bool front_layer_change = false;
            int phy_q1,phy_q2;
            int two_bit_max_decay;
            for (auto gate:comb) {               
                phy_q1 = gate->control;
                phy_q2 = gate->target;

                two_bit_max_decay = std::max(this->decay[phy_q1], this->decay[phy_q2]);
                this->decay[phy_q1] = two_bit_max_decay + 3*this->chip->coupling_map[phy_q1][phy_q2];
                this->decay[phy_q2] = this->decay[phy_q1];
                this->max_decay = std::max(this->decay[phy_q1], this->max_decay);
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
                    continue;
                }
    

                
                two_bit_max_decay = std::max(this->real_decay[phy_q1], this->real_decay[phy_q2]);
                this->real_decay[phy_q1] = two_bit_max_decay + 3*this->chip->coupling_map[phy_q1][phy_q2];
                this->real_decay[phy_q2] = this->real_decay[phy_q1];
                
                run_H_gate_before_op(gate->control, two_bit_max_decay);
                run_H_gate_before_op(gate->target, two_bit_max_decay);

                auto swap_gate = Gate(SWAP, 0, phy_q1, phy_q2);
                this->last_swap_pair[phy_q1] = phy_q2;
                this->last_swap_pair[phy_q2] = phy_q1;
                output.push_back(std::move(swap_gate));
                this->last_swap_gate[phy_q1].push_back(&output.back());
                this->last_swap_gate[phy_q2].push_back(&output.back());
                
            }
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
                    output.push_back(Gate(tmp->type, tmp->id, phy_qubit, -1));
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
            output.push_back(Gate(CX,node->cx_gate->id, phy_q1, phy_q2));
        }

        void update_front_layer() {

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
                    extended_layer.push_back(succ);
                    succ->has_done = true;
                }
            }

            for (auto node:extended_layer)
                node->has_done = false;        
        
            //同步decay 与real_decay;
            this->max_decay = 0;
            for(int i = 0; i<this->chip->qubit_number; i++) {
                decay[i] = real_decay[i];
                max_decay = std::max(max_decay, decay[i]);
            }

            // vector<int> succ_cnt = vector<int>(front_layer.size(), 0);
            // for(int i = 0; i<front_layer.size(); i++) {
            //     for (auto succ:front_layer[i]->successors) {
            //         if (succ->front_cnt == 1) succ_cnt[i] ++;
            //     }
            // }
            //排序front layer
            std::sort(this->front_layer.begin(), this->front_layer.end(), [this](Node* i, Node *j){

                int i_succ_cnt = 0;
                int j_succ_cnt = 0;
                for(auto succ:i->successors) if (succ->layer == 1) i_succ_cnt ++;
                for(auto succ:j->successors) if (succ->layer == 1) j_succ_cnt ++;
                return i_succ_cnt > j_succ_cnt;
                //return i->successors.size() > j->successors.size();
            });

        }

        

};


int main(int argc, char **argv) {
    if (argc != 7) {
        return 1;
    }
    char *chip_path = argv[1];
    char *circuit_path = argv[2];
    char *result_path = argv[3];
    DECAY_WEIGHT = std::stof(argv[4]);
    SWAP_PRIOR_TH = std::stof(argv[5]);
    NODE_PRIOR_TH = std::stof(argv[6]);

    Circuit* c = new Circuit(circuit_path);
    Quantum_chip *chip = new Quantum_chip("aaa", chip_path);
    vector<int> initial_layout(c->qubit_number);
    for (int i =0; i<c->qubit_number; i++)
        initial_layout[i] = i;
    Mapper mp(chip, c, initial_layout);
    clock_t start_time = clock();
    mp.route();
    clock_t end_time = clock();
    double run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
    std::cout<<run_time<<std::endl;
    mp.dump_result(result_path);
    delete c;
    delete chip;
    return 0;
}