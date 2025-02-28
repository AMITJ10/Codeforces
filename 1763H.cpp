# include <iostream>
# include <set>
# include <vector>
# include <algorithm>

using namespace std;

constexpr int kMaxNodes = 2e5 + 7;
constexpr long long kModulo = 998244353;
constexpr long long kInf = 1e18 + 7;
using NodeSet = pair<vector<long long>, long long>;

long long numNodes, queries, dfsTimer, logTable[kMaxNodes], depth[kMaxNodes], power[kMaxNodes];
long long dfsPos[kMaxNodes], reverseDfs[kMaxNodes], minTable[22][kMaxNodes], nodeValues[kMaxNodes];

vector<vector<long long>> adjacency(kMaxNodes);
set<NodeSet> activeNodes;

/**
 * @brief Finds the node with the minimum depth between two given nodes.
 *
 * @param nodeA The first node.
 * @param nodeB The second node.
 * @return long long The node with the minimum depth.
 */
long long getMinDepthNode(long long nodeA, long long nodeB) {
    return depth[nodeA] < depth[nodeB] ? nodeA : nodeB;
}

/**
 * @brief Finds the node with the maximum depth between two given nodes.
 *
 * @param nodeA The first node.
 * @param nodeB The second node.
 * @return long long The node with the maximum depth.
 */
long long getMaxDepthNode(long long nodeA, long long nodeB) {
    return depth[nodeA] > depth[nodeB] ? nodeA : nodeB;
}

/**
 * @brief Performs Depth-First Search (DFS) to populate depth and LCA-related structures.
 *
 * @param node The current node.
 * @param parent The parent node.
 */
void performDfs(long long node, long long parent) {
    reverseDfs[dfsPos[node] = ++dfsTimer] = node;
    minTable[0][dfsTimer] = parent;
    depth[node] = depth[parent] + 1;

    for (auto& child : adjacency[node]) {
        performDfs(child, node);
    }
}

/**
 * @brief Finds the Lowest Common Ancestor (LCA) of two nodes using RMQ technique.
 *
 * @param node1 The first node.
 * @param node2 The second node.
 * @return long long The LCA of the two nodes.
 */
long long findLca(long long node1, long long node2) {
    if (node1 == node2) return node1;

    long long pos1 = dfsPos[node1], pos2 = dfsPos[node2];
    if (pos1 > pos2) swap(pos1, pos2);

    long long logVal = logTable[pos2 - pos1++];
    return getMinDepthNode(minTable[logVal][pos1], minTable[logVal][pos2 - power[logVal] + 1]);
}

/**
 * @brief Computes the LCA-related value for an iterator in the active nodes set.
 *
 * @param it The iterator pointing to the current node in the active nodes set.
 * @return long long The computed LCA-related value.
 */
long long computeLcaResult(auto it) {
    auto nextIt = it, prevIt = it;
    long long result = 0;

    if ((nextIt = next(it))->first == it->first) {
        result = getMaxDepthNode(result, findLca(reverseDfs[it->second], reverseDfs[nextIt->second]));
    }
    if ((prevIt = prev(it))->first == it->first) {
        result = getMaxDepthNode(result, findLca(reverseDfs[it->second], reverseDfs[prevIt->second]));
    }
    return result;
}

/**
 * @brief Inserts a node into the active node set.
 *
 * @param node The node to be inserted.
 */
void insertNode(long long node) {
    if (!nodeValues[node]) return;

    vector<long long> temp;
    for (int i = 30; i >= 0; --i) {
        if ((nodeValues[node] - 1) >> i & 1) {
            temp.push_back(depth[node] + i);
        }
    }

    auto it = activeNodes.insert(NodeSet(temp, dfsPos[node])).first;
    while (it != activeNodes.end()) {
        long long res = computeLcaResult(it);
        if (res == 0) break;
        temp.push_back(depth[res]);
        it = activeNodes.insert(NodeSet(temp, dfsPos[res])).first;
    }
}

/**
 * @brief Removes a node from the active node set.
 *
 * @param node The node to be removed.
 */
void eraseNode(long long node) {
    if (!nodeValues[node]) return;

    vector<long long> temp;
    for (int i = 30; i >= 0; --i) {
        if ((nodeValues[node] - 1) >> i & 1) {
            temp.push_back(depth[node] + i);
        }
    }

    auto it = activeNodes.find(NodeSet(temp, dfsPos[node]));
    while (it != activeNodes.end()) {
        long long res = computeLcaResult(it);
        activeNodes.erase(it);
        if (res == 0) break;
        temp.push_back(depth[res]);
        it = activeNodes.find(NodeSet(temp, dfsPos[res]));
    }
}

/**
 * @brief Calculates the final result based on the active node set.
 *
 * @return long long The calculated final result.
 */
long long calculateFinalResult() {
    if (activeNodes.size() == 2) return 0;

    long long result = 1;
    for (auto i : (--(--activeNodes.end()))->first) {
        (result += power[i]) %= kModulo;
    }
    return result;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> numNodes;
    depth[numNodes + 1] = kInf;
    depth[0] = -1;

    activeNodes.insert(NodeSet(vector<long long>{ kInf }, 0));
    activeNodes.insert(NodeSet(vector<long long>{ -1 }, 0));

    for (long long i = 2; i <= numNodes; ++i) {
        long long parent;
        cin >> parent;
        adjacency[parent].push_back(i);
    }

    power[0] = 1;
    for (long long i = 1; i < kMaxNodes; ++i) {
        power[i] = power[i - 1] * 2 % kModulo;
    }

    logTable[0] = -1;
    for (long long i = 1; i < kMaxNodes; ++i) {
        logTable[i] = logTable[i >> 1] + 1;
    }

    performDfs(1, 0);

    for (long long i = 1; i <= 20; ++i) {
        for (long long j = 1; j + power[i] - 1 <= numNodes; ++j) {
            minTable[i][j] = getMinDepthNode(minTable[i - 1][j], minTable[i - 1][j + power[i - 1]]);
        }
    }

    for (long long i = 1; i <= numNodes; ++i) {
        cin >> nodeValues[i];
        insertNode(i);
    }

    cout << calculateFinalResult() << '\n';

    cin >> queries;
    while (queries--) {
        long long node, value;
        cin >> node >> value;
        eraseNode(node);
        nodeValues[node] = value;
        insertNode(node);
        cout << calculateFinalResult() << '\n';
    }

    return 0;
}
