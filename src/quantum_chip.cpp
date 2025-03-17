#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <string.h>
#include <cstdint>
#include <cmath>
#include "quantum_chip.h"

Quantum_chip::Quantum_chip(std::string name, char * file_path):
    name(name), gate_id(std::unordered_map<std::string, int>()), min_cx(INT32_MAX), max_cx(0),max_T1(0), max_T2(0), min_h(INT32_MAX){
    if (this->load(file_path) == 0) {
        this->calc_dist_m();
        this->neibors = std::vector<std::vector<int>>
                        (this->qubit_number,std::vector<int>());
        for (int i = 0; i<this->qubit_number; i++){
            for (int j = 0; j<this->qubit_number; j++){
                if (this->coupling_map[i][j] != 0)
                    this->neibors[i].push_back(j);
            }
        }
    }
}
// Quantum_chip::~Quantum_chip(){
// }


int Quantum_chip::load(const char * filepath){
    int fd = open(filepath, O_RDONLY);
    if (fd == -1) {
        std::cerr << "Could not open file." << std::endl;
        return 1;
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        std::cerr << "Could not get file size." << std::endl;
        close(fd);
        return 1;
    }

    size_t fileSize = sb.st_size;

    void* fileData = mmap(NULL, fileSize, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
    if (fileData == MAP_FAILED) {
        std::cerr << "Could not map file to memory." << std::endl;
        close(fd);
        return 0;
    }
    // 关闭文件描述符，数据仍然在内存中
    close(fd);

    // 读取文件数据
    char* data = static_cast<char*>(fileData);
    char* p = data;
    qubit_number = strtol(p, &p, 10);
    p++;
    coupling_map = std::vector<std::vector<int>>
                    (qubit_number, std::vector<int>(qubit_number, 0));
    
    double_gate_fidelity = std::vector<std::vector<int>>
                    (qubit_number, std::vector<int>(qubit_number, 0));
    qubit_T1 = std::vector<int>(qubit_number, 0);
    qubit_T2 = std::vector<int>(qubit_number, 0);

    int coupling_number = strtol(p, &p, 10);
    p++;
    //skip these data
    int row_cnt = 0;
    while (row_cnt<coupling_number) {
        while (*p != '\n') p++;
        p++;
        row_cnt++;
    }

    int T1,T2;
    for(int i = 0; i<this->qubit_number; i++) {
        T1 = strtol(p,&p,10);
        while (*p == ' ') p++;
        T2 = strtol(p,&p,10);
        while (*p != '\n') p++;
        p++;

        this->qubit_T1[i] = T1;
        this->qubit_T2[i] = T2;
        
    }
    

    char *z;
    //read single gate set;
    
    while (*p != '\n') {
        while (*p == ' ') p++;
        z = p;
        while (*p != ' ' && *p != '\n') p++;
        gate_set.push_back(std::string(z, p-z));
    }
    p++;
    int gate_set_number = gate_set.size();
    double_gate = gate_set[gate_set_number-1];

    single_gate_time = std::vector<std::vector<int>>
                        (qubit_number, std::vector<int>(gate_set_number-1, 0));
    single_gate_fidelity = std::vector<std::vector<int>>
                        (qubit_number, std::vector<int>(gate_set_number-1, 0));

    for (int i = 0; i<gate_set_number; i++) {
        gate_id[gate_set[i]] = i;
    }

    //std::string tmpgate;
    int bit1, bit2;
    float error;
    int gatetime;
    while (p<data+fileSize) {
        z = p;
        while (*p != ' ' && *p != '\n') p++;
        if (*p == '\n') {
            p++;
            continue;
        }
        *p = '\0';
        p++;
        while (*p == ' ') p++;
        //判断是否是单门还是双门
        if (strcmp(z, double_gate.c_str())) {
            //单门
            bit1 = strtol(p, &p, 10);
            while (*p == ' ') p++;
            gatetime = (int)(strtof(p, &p));
            while (*p == ' ') p++;
            error = strtof(p, &p);
            single_gate_time[bit1][gate_id[z]] = gatetime;
            single_gate_fidelity[bit1][gate_id[z]] = (int)(-std::log10(1-error)*1000000);
        } else {
            //双门
            bit1 = strtol(p, &p, 10);
            while (*p == ' ') p++;
            bit2 = strtol(p, &p, 10);
            while (*p == ' ') p++;
            gatetime = (int)(strtof(p, &p));
            while (*p == ' ') p++;

            error = strtof(p, &p);

            coupling_map[bit1][bit2] = gatetime;
            coupling_map[bit2][bit1] = gatetime;

            double_gate_fidelity[bit1][bit2] = (int)(-std::log10(1-error)*1000000);
            double_gate_fidelity[bit2][bit1] = double_gate_fidelity[bit1][bit2];
        }
        p++;
    }

    munmap(fileData, fileSize);
    return 0;
}

void Quantum_chip::calc_dist_m() {
    dist = std::vector<std::vector<int>>
            (qubit_number, std::vector<int>(qubit_number, INT32_MAX));
    steps_dist = std::vector<std::vector<int>>
            (qubit_number, std::vector<int>(qubit_number, INT32_MAX));
    fidelity_dist = std::vector<std::vector<int>>
            (qubit_number, std::vector<int>(qubit_number, INT32_MAX));
    // next_path = std::vector<std::vector<int>>
            //  (qubit_number, std::vector<int>(qubit_number, 0));
    
    //初始化
    for (int i = 0; i < qubit_number; ++i) {
        dist[i][i] = 0;
        steps_dist[i][i] = 0;
        fidelity_dist[i][i] = 0;
        for (int j = 0; j < qubit_number; ++j)
            if (coupling_map[i][j] != 0) {
                dist[i][j] = coupling_map[i][j];
                steps_dist[i][j] = 200;
                fidelity_dist[i][j] = double_gate_fidelity[i][j];
                //if (i!=j) next_path[i][j] = j;
            }
    }


    //求解带权最短距离
    for (int k = 0; k < qubit_number; ++k) {
        // i is the source vertex
        for (int i = 0; i < qubit_number; ++i) {
            // j is the destination vertex
            for (int j = 0; j < qubit_number; ++j) {
                if (dist[i][k] != INT32_MAX && dist[k][j] != INT32_MAX
                    && dist[i][k] + dist[k][j] < dist[i][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    //next_path[i][j] = next_path[i][k];
                }
            }
        }
    }

    //求解具体路径
    // for (int i = 0; i < qubit_number; ++i) {
    //     for (int j = 0; j < qubit_number; ++j) {
    //         int v = i;
    //         while (v != j) {
    //             shortest_path[i][j].push_back(v);
    //             v = next_path[v][j];
    //         }
    //         shortest_path[i][j].push_back(j);
    //     }
    // }
    

    //求解保真度最短距离
    for (int k = 0; k < qubit_number; ++k) {
        // i is the source vertex
        for (int i = 0; i < qubit_number; ++i) {
            // j is the destination vertex
            for (int j = 0; j < qubit_number; ++j) {
                if (fidelity_dist[i][k] != INT32_MAX && fidelity_dist[k][j] != INT32_MAX
                    && fidelity_dist[i][k] + fidelity_dist[k][j] < fidelity_dist[i][j])
                {
                    fidelity_dist[i][j] = fidelity_dist[i][k] + fidelity_dist[k][j];
                }
            }
        }
    }



    //求解无权最短距离
    for (int k = 0; k < qubit_number; ++k) {
        // i is the source vertex
        for (int i = 0; i < qubit_number; ++i) {
            // j is the destination vertex
            for (int j = 0; j < qubit_number; ++j) {
                if (steps_dist[i][k] != INT32_MAX && steps_dist[k][j] != INT32_MAX
                    && steps_dist[i][k] + steps_dist[k][j] < steps_dist[i][j])
                {
                    steps_dist[i][j] = steps_dist[i][k] + steps_dist[k][j];
                }
            }
        }
    }

    //average_h_time_in_path = std::vector<std::vector<int>>
    //    (qubit_number, std::vector<int>(qubit_number, 0));
    // int cnt, total, v;
    //得出一些最值
    for (int i = 0; i<qubit_number; i++){
        for (int j = i; j<qubit_number; j++){
            // v = i;
            // cnt = 1;
            // //total = single_gate_time[i][0];
            // while (v != j) {
            //     //v = next_path[v][j];
            //     //cnt +=1;
            //     //total += single_gate_time[v][0];
            // }
            // //average_h_time_in_path[i][j] = total/cnt;
            // //average_h_time_in_path[j][i] = total/cnt;
            if (this->coupling_map[i][j] == 0) continue;
            this->min_cx = std::min(this->min_cx, this->coupling_map[i][j]);
            this->max_cx = std::max(this->max_cx, this->coupling_map[i][j]);
            
        }
        this->max_T1 = std::max(this->max_T1, this->qubit_T1[i]);
        this->max_T2 = std::max(this->max_T2, this->qubit_T2[i]);
        //this->min_h = std::min(this->min_h, this->single_gate_time[i][0]);
    }
    
}