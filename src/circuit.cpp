#include "circuit.h"
#include <fstream>
#include <regex>
#include <iostream>

//******** Node */

void Node::add_successors(Node *n) {
    this->successors.push_back(n);
    if (this->layer >= n->layer) n->layer = this->layer + 1;
    n->front_cnt += 1;
}


//******* Circuit */

Circuit::Circuit() {
    qubit_number = 0;
}

Circuit::Circuit(char *file_path) {
    qubit_number = 0;
    this->load(file_path);
}

bool Circuit::load(char* file_path) {
    std::ifstream file(file_path); // 打开文件
    if (!file.is_open()) {
        std::cout << "Failed to open file." << std::endl;
        return false;
    }

    std::string line;
    std::smatch matches;
    while (std::getline(file, line)) {
        std::regex qreg_pattern(R"(\wreg \w\[(\d+)\];)");
        if (std::regex_match(line, matches, qreg_pattern)) {
            if (line[0] == 'q')
                this->qubit_number = std::stoi(matches[1].str());
            else {
                std::getline(file, line);
                break;
            }
        } else if (this->qubit_number > 0) break;
    }

    //this->qubit_gate_list = vector<vector<Gate *>>(qubit_number, vector<Gate *>());

    std::regex pattern(R"(.+? q\[(\d+)\](?:,q\[(\d+)\])?;)");
    Gate* gate;
    Node* node;
    Node* qbit_last_node[this->qubit_number] = {nullptr};
    logical_bit_node_list = vector<vector<int>>(this->qubit_number, vector<int>());
    vector<vector<Gate *>> tmp_h_gate_list = vector<vector<Gate *>>(this->qubit_number, vector<Gate *>());
    do { // 逐行读取文件
        if (std::regex_match(line, matches, pattern)) { // 匹配每一行
            //std::cout << "Match found: " << matches[0] << std::endl;
            
            // 将捕获的数字字符串转换为整数
            int q1 = std::stoi(matches[1].str());

            // 检查并处理第二个捕获组
            if (matches[2].matched) {
                int q2 = std::stoi(matches[2].str());
                gate = new Gate(CX, this->gate_list.size(), q1, q2);
                node = new Node(this->node_list.size());
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
                this->node_list.push_back(node);
                logical_bit_node_list[q1].push_back(node->id);
                logical_bit_node_list[q2].push_back(node->id);
                //std::cout << "CX " << q1 << q2 << std::endl;
                
            } else {
                //std::cout<< "H " << q1 << std::endl;
                if (line[0] == 't')
                    gate = new Gate(T, this->gate_list.size(), q1, -1);
                else
                    gate = new Gate(H, this->gate_list.size(), q1, -1);
                tmp_h_gate_list[q1].push_back(gate);
            }
            this->gate_list.push_back(gate);
            //this->qubit_gate_list[gate->control].push_back(gate);
            //if (gate->target != -1) this->qubit_gate_list[gate->target].push_back(gate);
        }
    } while (std::getline(file, line));
    this->left_h_gate_each_bit.swap(tmp_h_gate_list);
    file.close(); // 关闭文件
    return true;
}

Circuit& Circuit::operator=(const Circuit & other) {
    if (this != &other) {
        logical_bit_node_list = other.logical_bit_node_list;
        for (size_t i = 0; i<node_list.size(); i++) {
            node_list[i]->has_done = other.node_list[i]->has_done;
            node_list[i]->front_cnt = other.node_list[i]->front_cnt;
            node_list[i]->layer = other.node_list[i]->layer;
            node_list[i]->control_h_to_do_index = other.node_list[i]->control_h_to_do_index;
            node_list[i]->target_h_to_do_index = other.node_list[i]->target_h_to_do_index;
        }
    }
    return *this;
}

Circuit::~Circuit(){
    for(Gate* gate:this->gate_list) delete gate;
}
