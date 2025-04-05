#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <iostream>
#include "quantum_chip.h"
#include <iostream>
#include <vector>
#include <algorithm> // for std::shuffle
#include <random>    // for std::mt19937
#include <chrono>    // for std::chrono::system_clock

static int W = 1;

void check(const std::vector<std::vector<int>>& cmapMatrix, 
    const std::vector<int>& originPbit2Lbit, 
    const std::vector<int>& targetPbit2Lbit,
    const std::vector<std::vector<int>> swapSequence) {
    
    auto origin = originPbit2Lbit;
    auto target = targetPbit2Lbit;
    for (auto swap:swapSequence) {
        int swapPhyInd1 = swap[0];
        int swapPhyInd2 = swap[1];
        if (cmapMatrix[swapPhyInd1][swapPhyInd2] == 0) {
            std::cout << "error" << std::endl;
            return;
        }
        int temp = origin[swapPhyInd1];
        origin[swapPhyInd1] = origin[swapPhyInd2];
        origin[swapPhyInd2] = temp;
    }
    for (int i = 0; i < origin.size(); i++) {
        if (origin[i] != target[i]) {
            std::cout << "error" << std::endl;
            return;
        }
    }
}

// 图结构定义
struct Graph {
    // 将邻接表改为存储整数边权重
    std::unordered_map<int, std::unordered_map<int, int>> adjacency;
    std::unordered_map<int, std::unordered_map<std::string, int>> nodes;
    
    void add_node(int node) {
        if (adjacency.find(node) == adjacency.end()) {
            adjacency[node] = std::unordered_map<int, int>();
        }
        if (nodes.find(node) == nodes.end()) {
            nodes[node] = std::unordered_map<std::string, int>();
        }
    }
    
    void add_nodes_from(const std::vector<int>& node_list) {
        for (int node : node_list) {
            add_node(node);
        }
    }
    
    void add_edge(int u, int v, int weight = 1) {
        add_node(u);
        add_node(v);
        adjacency[u][v] = weight;
        adjacency[v][u] = weight;
    }
    
    // 从邻接矩阵添加边
    void add_edges_from_matrix(const std::vector<std::vector<int>>& adjacency_matrix) {
        for (int i = 0; i < adjacency_matrix.size(); ++i) {
            for (int j = i + 1; j < adjacency_matrix[i].size(); ++j) {
                if (adjacency_matrix[i][j] != 0) {
                    add_edge(i, j, 1);
                }
            }
        }
    }
};

struct DiGraph {
    std::unordered_map<int, std::unordered_map<int, int>> adjacency;
    std::unordered_map<int, std::unordered_map<std::string, int>> nodes;
    
    void add_node(int node) {
        if (adjacency.find(node) == adjacency.end()) {
            adjacency[node] = std::unordered_map<int, int>();
        }
        if (nodes.find(node) == nodes.end()) {
            nodes[node] = std::unordered_map<std::string, int>();
        }
    }
    
    void add_nodes_from(const std::vector<int>& node_list) {
        for (int node : node_list) {
            add_node(node);
        }
    }
    
    void add_edge(int u, int v, int weight = 1) {
        add_node(u);
        add_node(v);
        adjacency[u][v] = weight;
    }
};

// 计算所有节点对之间的最短路径长度
std::unordered_map<int, std::unordered_map<int, int>> all_pairs_shortest_path_length(const Graph& graph) {
    std::unordered_map<int, std::unordered_map<int, int>> distance;
    
    for (const auto& node_pair : graph.adjacency) {
        int source = node_pair.first;
        std::unordered_map<int, int> dist;
        std::unordered_set<int> visited;
        
        // 初始化距离
        for (const auto& adj_pair : graph.adjacency) {
            dist[adj_pair.first] = INT32_MAX;
        }
        dist[source] = 0;
        
        while (visited.size() < graph.adjacency.size()) {
            // 找到未访问节点中距离最小的
            int current = -1;
            int min_dist = INT32_MAX;
            for (const auto& d : dist) {
                if (visited.find(d.first) == visited.end() && d.second < min_dist) {
                    current = d.first;
                    min_dist = d.second;
                }
            }
            
            if (current == -1) break;
            visited.insert(current);
            
            // 更新邻居距离
            for (const auto& neighbor : graph.adjacency.at(current)) {
                int next = neighbor.first;
                int weight = neighbor.second;
                if (dist[next] > dist[current] + weight) {
                    dist[next] = dist[current] + weight;
                }
            }
        }
        
        distance[source] = dist;
    }
    
    return distance;
}

// 构建有向图
DiGraph ConstructDirectedGraph(
    const std::unordered_map<int, std::unordered_map<int, int>>& pDistance, 
    const Graph& undirectedGraph, 
    const std::vector<int>& targetPbit2Lbit, 
    const std::unordered_map<int, int>& targetLbit2Pbit, 
    const std::unordered_map<int, std::unordered_map<int, int>>& adjacency){
    
    DiGraph dGraph;
    std::vector<int> nodes;
    

    for (const auto& node_pair : undirectedGraph.nodes) {
        nodes.push_back(node_pair.first);
    }
    

    if (W == 1) {
        // 对nodes进行排序
        std::sort(nodes.begin(), nodes.end(), [&](int a, int b) {
            return pDistance.at(a).at(targetLbit2Pbit.at(undirectedGraph.nodes.at(a).at("lbit")))
            > pDistance.at(b).at(targetLbit2Pbit.at(undirectedGraph.nodes.at(b).at("lbit")));
        });
    }


    dGraph.add_nodes_from(nodes);
    
    for (int node : nodes) {
        if (undirectedGraph.nodes.at(node).at("lbit") != targetPbit2Lbit[node]) {
            for (const auto& adj_pair : adjacency.at(node)) {
                int adjNode = adj_pair.first;
                int targetPbitOfNode = targetLbit2Pbit.at(undirectedGraph.nodes.at(node).at("lbit"));
                if (pDistance.at(adjNode).at(targetPbitOfNode) < pDistance.at(node).at(targetPbitOfNode)) {
                    dGraph.add_edge(node, adjNode);
                }
            }
        }
    }
    
    return dGraph;
}

std::pair<int, int> get_next_swap(
    const std::unordered_map<int, std::unordered_map<int, int>>& pDistance, 
    const Graph& undirectedGraph, 
    const std::vector<int>& targetPbit2Lbit, 
    const std::unordered_map<int, int>& targetLbit2Pbit, 
    const std::unordered_map<int, std::unordered_map<int, int>>& adjacency) {
    
    for (const auto& node_pair : undirectedGraph.nodes) {
        int node = node_pair.first;
        //if (undirectedGraph.nodes.at(node).at("lbit") != targetPbit2Lbit[node]) {
            for (const auto& adj_pair : adjacency.at(node)) {
                int adjNode = adj_pair.first;
                int targetPbitOfNode = targetLbit2Pbit.at(undirectedGraph.nodes.at(node).at("lbit"));
                int targetPbitOfAdjNode = targetLbit2Pbit.at(undirectedGraph.nodes.at(adjNode).at("lbit"));

                int deltaofNodeLbit = pDistance.at(adjNode).at(targetPbitOfNode) - pDistance.at(node).at(targetPbitOfNode) < 0;
                int deltaofadjNodeLbit = pDistance.at(node).at(targetPbitOfAdjNode) - pDistance.at(adjNode).at(targetPbitOfAdjNode) < 0;

                if ((deltaofNodeLbit && deltaofadjNodeLbit) ||
                    (deltaofNodeLbit && adjNode == targetPbitOfAdjNode) ||
                    (node == targetPbitOfNode && deltaofadjNodeLbit)) {
                    return std::make_pair(node, adjNode);
                }
            }
        //}
    }
    return std::make_pair(-1, -1); // No swap found
}

// TS4 算法实现，使用邻接矩阵表示耦合图
std::vector<std::vector<int>> TS4(
    const std::vector<std::vector<int>>& cmapMatrix, 
    const std::vector<int>& originPbit2Lbit, 
    const std::vector<int>& targetPbit2Lbit,
    double *run_time) {
    
    // 构建物理图
    Graph physicalGraph;
    
    // 添加节点
    std::vector<int> nodeList;
    for (int i = 0; i < cmapMatrix.size(); ++i) {
        nodeList.push_back(i);
    }
    physicalGraph.add_nodes_from(nodeList);
    
    // 从邻接矩阵添加边
    physicalGraph.add_edges_from_matrix(cmapMatrix);
    
    // 计算所有节点对之间的最短路径长度
    auto pDistance = all_pairs_shortest_path_length(physicalGraph);
    
    // 创建目标逻辑比特到物理比特的映射
    std::unordered_map<int, int> targetLbit2Pbit;
    int notInPositionCount = 0;
    
    // 初始化节点属性和映射
    for (int pIndex = 0; pIndex < targetPbit2Lbit.size(); pIndex++) {
        physicalGraph.nodes[pIndex]["lbit"] = originPbit2Lbit[pIndex];
        targetLbit2Pbit[targetPbit2Lbit[pIndex]] = pIndex;
        if (originPbit2Lbit[pIndex] != targetPbit2Lbit[pIndex]) {
            notInPositionCount++;
        }
    }
    
    std::vector<std::vector<int>> swapSequence;
    clock_t start_time = clock();

    while (notInPositionCount > 0) {
        DiGraph assistDiGraph = ConstructDirectedGraph(
            pDistance, physicalGraph, targetPbit2Lbit, targetLbit2Pbit, physicalGraph.adjacency);
        
        std::vector<int> cycleNodes;
        
        for (const auto &node_info : assistDiGraph.nodes) {
            int node = node_info.first;  // 改为遍历排序后的nodes
            if (assistDiGraph.adjacency[node].empty()) {
                continue;
            } else {
                cycleNodes.push_back(node);
                int nodePointer = node;
                
                while (!assistDiGraph.adjacency[nodePointer].empty()) {
                    std::vector<int> possibleNextNodeList;
                    for (const auto& next_pair : assistDiGraph.adjacency[nodePointer]) {
                        possibleNextNodeList.push_back(next_pair.first);
                    }
                    
                    nodePointer = -1;
                    for (int candidate : possibleNextNodeList) {
                        if (nodePointer == -1 || 
                            physicalGraph.nodes[nodePointer]["lbit"] == targetPbit2Lbit[nodePointer]) {
                            nodePointer = candidate;
                        }
                    }
                    
                    if (std::find(cycleNodes.begin(), cycleNodes.end(), nodePointer) != cycleNodes.end()) {
                        break;
                    } else {
                        cycleNodes.push_back(nodePointer);
                    }
                }
                
                if (assistDiGraph.adjacency[nodePointer].empty()) {
                    int swapPhyInd1 = cycleNodes.back();
                    int swapPhyInd2 = cycleNodes[cycleNodes.size() - 2];
                    swapSequence.push_back({swapPhyInd1, swapPhyInd2});
                    
                    int temp = physicalGraph.nodes[swapPhyInd1]["lbit"];
                    physicalGraph.nodes[swapPhyInd1]["lbit"] = physicalGraph.nodes[swapPhyInd2]["lbit"];
                    physicalGraph.nodes[swapPhyInd2]["lbit"] = temp;
                    notInPositionCount += 1;
                } else {
                    int swapInd = -1;
                    while (cycleNodes[cycleNodes.size() + swapInd] != nodePointer) {
                        int swapPhyInd1 = cycleNodes[cycleNodes.size() + swapInd];
                        int swapPhyInd2 = cycleNodes[cycleNodes.size() + swapInd - 1];
                        swapSequence.push_back({swapPhyInd1, swapPhyInd2});
                        
                        int temp = physicalGraph.nodes[swapPhyInd1]["lbit"];
                        physicalGraph.nodes[swapPhyInd1]["lbit"] = physicalGraph.nodes[swapPhyInd2]["lbit"];
                        physicalGraph.nodes[swapPhyInd2]["lbit"] = temp;
                        
                        if (physicalGraph.nodes[swapPhyInd1]["lbit"] == targetPbit2Lbit[swapPhyInd1]) {
                            notInPositionCount--;
                        }
                        swapInd--;
                    }
                    
                    if (physicalGraph.nodes[cycleNodes[cycleNodes.size() + swapInd]]["lbit"] == 
                        targetPbit2Lbit[cycleNodes[cycleNodes.size() + swapInd]]) {
                        notInPositionCount--;
                    }
                }
                break;
            }
        }
    }
    clock_t end_time = clock();
    *run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
    //check(cmapMatrix, originPbit2Lbit, targetPbit2Lbit, swapSequence);
    return swapSequence;
}

std::vector<std::vector<int>> TS4_new(
    const std::vector<std::vector<int>>& cmapMatrix, 
    const std::vector<int>& originPbit2Lbit, 
    const std::vector<int>& targetPbit2Lbit,
    double *run_time) {
    
    // 构建物理图
    Graph physicalGraph;
    
    // 添加节点
    std::vector<int> nodeList;
    for (int i = 0; i < cmapMatrix.size(); ++i) {
        nodeList.push_back(i);
    }
    physicalGraph.add_nodes_from(nodeList);
    
    // 从邻接矩阵添加边
    physicalGraph.add_edges_from_matrix(cmapMatrix);
    
    // 计算所有节点对之间的最短路径长度
    auto pDistance = all_pairs_shortest_path_length(physicalGraph);
    
    // 创建目标逻辑比特到物理比特的映射
    std::unordered_map<int, int> targetLbit2Pbit;
    int notInPositionCount = 0;
    
    // 初始化节点属性和映射
    for (int pIndex = 0; pIndex < targetPbit2Lbit.size(); pIndex++) {
        physicalGraph.nodes[pIndex]["lbit"] = originPbit2Lbit[pIndex];
        targetLbit2Pbit[targetPbit2Lbit[pIndex]] = pIndex;
        if (originPbit2Lbit[pIndex] != targetPbit2Lbit[pIndex]) {
            notInPositionCount++;
        }
    }
    
    std::vector<std::vector<int>> swapSequence;
    clock_t start_time = clock();

    while (notInPositionCount > 0) {        
        auto swap = get_next_swap(
            pDistance, physicalGraph, targetPbit2Lbit, targetLbit2Pbit, physicalGraph.adjacency);
        int swapPhyInd1 = swap.first;
        int swapPhyInd2 = swap.second;
        
        if (swapPhyInd1 != -1 && swapPhyInd2 != -1) {
            //std::cout << swapPhyInd1 << "\t\t" << swapPhyInd2 << std::endl;
            // 执行一次交换
            swapSequence.push_back({swapPhyInd1, swapPhyInd2});
            
            // 交换逻辑比特
            int temp = physicalGraph.nodes[swapPhyInd1]["lbit"];
            physicalGraph.nodes[swapPhyInd1]["lbit"] = physicalGraph.nodes[swapPhyInd2]["lbit"];
            physicalGraph.nodes[swapPhyInd2]["lbit"] = temp;
            
            // 更新计数
            if (physicalGraph.nodes[swapPhyInd1]["lbit"] == targetPbit2Lbit[swapPhyInd1]) {
                notInPositionCount--;
            }
            if (physicalGraph.nodes[swapPhyInd2]["lbit"] == targetPbit2Lbit[swapPhyInd2]) {
                notInPositionCount--;
            }
            //交换前就在正确位置，交换后被改变
            if (physicalGraph.nodes[swapPhyInd2]["lbit"] == targetPbit2Lbit[swapPhyInd1]) {
                notInPositionCount++;
            }
            if (physicalGraph.nodes[swapPhyInd1]["lbit"] == targetPbit2Lbit[swapPhyInd2]) {
                notInPositionCount++;
            }

        } else {
            DiGraph assistDiGraph = ConstructDirectedGraph(
                pDistance, physicalGraph, targetPbit2Lbit, targetLbit2Pbit, physicalGraph.adjacency);
            
            std::vector<int> cycleNodes;
            
            for (const auto &node_info : assistDiGraph.nodes) {
                int node = node_info.first;  // 改为遍历排序后的nodes
                if (assistDiGraph.adjacency[node].empty()) {
                    continue;
                } else {
                    cycleNodes.push_back(node);
                    int nodePointer = node;
                    
                    while (!assistDiGraph.adjacency[nodePointer].empty()) {
                        std::vector<int> possibleNextNodeList;
                        for (const auto& next_pair : assistDiGraph.adjacency[nodePointer]) {
                            possibleNextNodeList.push_back(next_pair.first);
                        }
                        
                        nodePointer = -1;
                        for (int candidate : possibleNextNodeList) {
                            if (nodePointer == -1 || 
                                physicalGraph.nodes[nodePointer]["lbit"] == targetPbit2Lbit[nodePointer]) {
                                nodePointer = candidate;
                            }
                        }
                        
                        if (std::find(cycleNodes.begin(), cycleNodes.end(), nodePointer) != cycleNodes.end()) {
                            break;
                        } else {
                            cycleNodes.push_back(nodePointer);
                        }
                    }
                    
                    if (assistDiGraph.adjacency[nodePointer].empty()) {
                        int swapPhyInd1 = cycleNodes.back();
                        int swapPhyInd2 = cycleNodes[cycleNodes.size() - 2];
                        swapSequence.push_back({swapPhyInd1, swapPhyInd2});
                        
                        int temp = physicalGraph.nodes[swapPhyInd1]["lbit"];
                        physicalGraph.nodes[swapPhyInd1]["lbit"] = physicalGraph.nodes[swapPhyInd2]["lbit"];
                        physicalGraph.nodes[swapPhyInd2]["lbit"] = temp;
                        notInPositionCount += 1;
                    } else {
                        int swapInd = -1;
                        while (cycleNodes[cycleNodes.size() + swapInd] != nodePointer) {
                            int swapPhyInd1 = cycleNodes[cycleNodes.size() + swapInd];
                            int swapPhyInd2 = cycleNodes[cycleNodes.size() + swapInd - 1];
                            swapSequence.push_back({swapPhyInd1, swapPhyInd2});
                            
                            int temp = physicalGraph.nodes[swapPhyInd1]["lbit"];
                            physicalGraph.nodes[swapPhyInd1]["lbit"] = physicalGraph.nodes[swapPhyInd2]["lbit"];
                            physicalGraph.nodes[swapPhyInd2]["lbit"] = temp;
                            
                            if (physicalGraph.nodes[swapPhyInd1]["lbit"] == targetPbit2Lbit[swapPhyInd1]) {
                                notInPositionCount--;
                            }
                            swapInd--;
                        }
                        
                        if (physicalGraph.nodes[cycleNodes[cycleNodes.size() + swapInd]]["lbit"] == 
                            targetPbit2Lbit[cycleNodes[cycleNodes.size() + swapInd]]) {
                            notInPositionCount--;
                        }
                    }
                    break;
                }
            }
        }
    }
    clock_t end_time = clock();
    *run_time = static_cast<double>(end_time-start_time)/ CLOCKS_PER_SEC;
    //check(cmapMatrix, originPbit2Lbit, targetPbit2Lbit, swapSequence);
    return swapSequence;
}
