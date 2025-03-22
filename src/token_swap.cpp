#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <iostream>

// 图结构定义
struct Graph {
    std::unordered_map<int, std::unordered_set<int>> adjacency;
    std::unordered_map<int, std::unordered_map<std::string, int>> nodes;
    
    void add_node(int node) {
        if (adjacency.find(node) == adjacency.end()) {
            adjacency[node] = std::unordered_set<int>();
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
    
    void add_edge(int u, int v) {
        add_node(u);
        add_node(v);
        adjacency[u].insert(v);
        adjacency[v].insert(u);
    }
    
    // 从邻接矩阵添加边
    void add_edges_from_matrix(const std::vector<std::vector<int>>& adjacency_matrix) {
        for (int i = 0; i < adjacency_matrix.size(); ++i) {
            for (int j = i + 1; j < adjacency_matrix[i].size(); ++j) {
                if (adjacency_matrix[i][j] != 0) {
                    add_edge(i, j);
                }
            }
        }
    }
};

struct DiGraph {
    std::unordered_map<int, std::unordered_set<int>> adjacency;
    std::unordered_map<int, std::unordered_map<std::string, int>> nodes;
    
    void add_node(int node) {
        if (adjacency.find(node) == adjacency.end()) {
            adjacency[node] = std::unordered_set<int>();
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
    
    void add_edge(int u, int v) {
        add_node(u);
        add_node(v);
        adjacency[u].insert(v);
    }
};

// 计算所有节点对之间的最短路径长度
std::unordered_map<int, std::unordered_map<int, int>> all_pairs_shortest_path_length(const Graph& graph) {
    std::unordered_map<int, std::unordered_map<int, int>> distance;
    
    for (const auto& node_pair : graph.adjacency) {
        int source = node_pair.first;
        std::unordered_map<int, int> dist;
        std::queue<int> q;
        
        dist[source] = 0;
        q.push(source);
        
        while (!q.empty()) {
            int current = q.front();
            q.pop();
            
            for (int neighbor : graph.adjacency.at(current)) {
                if (dist.find(neighbor) == dist.end()) {
                    dist[neighbor] = dist[current] + 1;
                    q.push(neighbor);
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
    const std::unordered_map<int, std::unordered_set<int>>& adjacency) {
    
    DiGraph dGraph;
    std::vector<int> nodes;
    
    for (const auto& node_pair : undirectedGraph.nodes) {
        nodes.push_back(node_pair.first);
    }
    
    dGraph.add_nodes_from(nodes);
    
    for (int node : nodes) {
        if (undirectedGraph.nodes.at(node).at("lbit") != targetPbit2Lbit[node]) {
            for (int adjNode : adjacency.at(node)) {
                int targetPbitOfNode = targetLbit2Pbit.at(undirectedGraph.nodes.at(node).at("lbit"));
                if (pDistance.at(adjNode).at(targetPbitOfNode) < pDistance.at(node).at(targetPbitOfNode)) {
                    dGraph.add_edge(node, adjNode);
                }
            }
        }
    }
    
    return dGraph;
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
        
        for (const auto& node_pair : assistDiGraph.adjacency) {
            int node = node_pair.first;
            if (node_pair.second.empty()) {
                continue;
            } else {
                cycleNodes.push_back(node);
                int nodePointer = node;
                
                while (!assistDiGraph.adjacency[nodePointer].empty()) {
                    std::vector<int> possibleNextNodeList;
                    for (int next : assistDiGraph.adjacency[nodePointer]) {
                        possibleNextNodeList.push_back(next);
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
    return swapSequence;
}

// 测试函数
// void test_TS4() {
//     // 定义邻接矩阵表示的耦合图
//     // 这里我们创建一个10x10的邻接矩阵，对应之前的耦合图
//     std::vector<std::vector<int>> cmapMatrix = {
//         {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
//         {1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
//         {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
//         {0, 0, 1, 0, 1, 0, 0, 0, 1, 0},
//         {0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
//         {0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
//         {0, 0, 0, 0, 0, 1, 0, 1, 0, 0},
//         {0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
//         {0, 0, 0, 1, 0, 0, 0, 1, 0, 1},
//         {0, 0, 0, 0, 1, 0, 0, 0, 1, 0}
//     };
    
//     std::vector<int> originLayout = {2, 3, 5, 7, 4, 9, 8, 0, 6, 1};
//     std::vector<int> targetLayout = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    
//     // 调用TS4算法
//     std::vector<std::vector<int>> swapSequence = TS4(cmapMatrix, originLayout, targetLayout);
    
//     // 打印结果
//     std::cout << "交换序列:" << std::endl;
//     for (const auto& swap : swapSequence) {
//         std::cout << "交换物理比特 " << swap[0] << " 和 " << swap[1] << std::endl;
//     }
    
//     // 验证结果
//     std::cout << "交换总数: " << swapSequence.size() << std::endl;
// }

// int main() {
//     test_TS4();
//     return 0;
// }