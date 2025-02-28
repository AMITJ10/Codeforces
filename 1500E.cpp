#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
using namespace std;

// Use `using` instead of `typedef` for better readability
using LongType = long long;

// Constants for array and segment tree sizes
constexpr int kMaxArraySize = 2e5 + 5;
constexpr int kMaxTreeSize = 4e5 + 5;

// Global variables
int array_size, query_count, unique_value_count;

/**
 * Represents a query with an operation type, an ID, and a value.
 */
struct Query {
    int operation_type;
    int id;
    LongType value;
} queries[kMaxArraySize];

vector<LongType> unique_values;
LongType array_values[kMaxArraySize], mapped_ids[kMaxArraySize];

/**
 * Maps a value to its unique ID in the sorted list of unique values.
 * @param value The value to map.
 * @return The unique ID of the value.
 */
int GetMappedId(LongType value) {
    return lower_bound(unique_values.begin(), unique_values.end(), value) - unique_values.begin() + 1;
}

/**
 * Represents a node in the segment tree.
 */
struct SegmentTreeNode {
    LongType sum;
    LongType count;
    LongType left_sum;
    LongType right_sum;
} segment_tree[kMaxTreeSize << 2];

/**
 * Updates the segment tree node by aggregating values from its children.
 * @param node The current node to update.
 */
void PushUp(int node) {
    segment_tree[node].count = segment_tree[node << 1].count + segment_tree[node << 1 | 1].count;
    segment_tree[node].sum = segment_tree[node << 1].sum + segment_tree[node << 1 | 1].sum;
    segment_tree[node].left_sum = segment_tree[node << 1].left_sum + segment_tree[node << 1 | 1].count * segment_tree[node << 1].sum + segment_tree[node << 1 | 1].left_sum;
    segment_tree[node].right_sum = segment_tree[node << 1 | 1].right_sum + segment_tree[node << 1].count * segment_tree[node << 1 | 1].sum + segment_tree[node << 1].right_sum;
}

/**
 * Updates the segment tree with a new value at a specific index.
 * @param node The current node.
 * @param left The left boundary of the current segment.
 * @param right The right boundary of the current segment.
 * @param index The index to update.
 * @param value The value to add or subtract.
 */
void UpdateSegmentTree(int node, int left, int right, int index, LongType value) {
    if (left == right) {
        segment_tree[node].count += (value < 0 ? -1 : 1);
        segment_tree[node].sum += value;
        segment_tree[node].left_sum += value;
        segment_tree[node].right_sum += value;
        return;
    }
    int mid = (left + right) >> 1;
    if (index <= mid)
        UpdateSegmentTree(node << 1, left, mid, index, value);
    else
        UpdateSegmentTree(node << 1 | 1, mid + 1, right, index, value);
    PushUp(node);
}

/**
 * Queries the sum of the leftmost `k` elements in the segment tree.
 * @param node The current node.
 * @param k The number of elements to query.
 * @return The sum of the leftmost `k` elements.
 */
LongType QueryLeft(int node, int k) {
    if (segment_tree[node].count <= k)
        return segment_tree[node].sum;
    LongType result = QueryLeft(node << 1, k);
    if (k > segment_tree[node << 1].count)
        result += QueryLeft(node << 1 | 1, k - segment_tree[node << 1].count);
    return result;
}

/**
 * Queries the sum of the rightmost `k` elements in the segment tree.
 * @param node The current node.
 * @param k The number of elements to query.
 * @return The sum of the rightmost `k` elements.
 */
LongType QueryRight(int node, int k) {
    if (segment_tree[node].count <= k)
        return segment_tree[node].sum;
    LongType result = QueryRight(node << 1 | 1, k);
    if (k > segment_tree[node << 1 | 1].count)
        result += QueryRight(node << 1, k - segment_tree[node << 1 | 1].count);
    return result;
}

/**
 * Queries the sum of the left sums for the leftmost `k` elements.
 * @param node The current node.
 * @param k The number of elements to query.
 * @return The sum of the left sums.
 */
LongType QuerySumLeft(int node, int k) {
    if (segment_tree[node].count <= k)
        return segment_tree[node].left_sum;
    LongType result = QuerySumLeft(node << 1, k);
    if (k > segment_tree[node << 1].count) {
        result += (k - segment_tree[node << 1].count) * segment_tree[node << 1].sum;
        result += QuerySumLeft(node << 1 | 1, k - segment_tree[node << 1].count);
    }
    return result;
}

/**
 * Queries the sum of the right sums for the rightmost `k` elements.
 * @param node The current node.
 * @param k The number of elements to query.
 * @return The sum of the right sums.
 */
LongType QuerySumRight(int node, int k) {
    if (segment_tree[node].count <= k)
        return segment_tree[node].right_sum;
    LongType result = QuerySumRight(node << 1 | 1, k);
    if (k > segment_tree[node << 1 | 1].count) {
        result += (k - segment_tree[node << 1 | 1].count) * segment_tree[node << 1 | 1].sum;
        result += QuerySumRight(node << 1, k - segment_tree[node << 1 | 1].count);
    }
    return result;
}

/**
 * Calculates the difference between the sum of the rightmost `x` elements
 * and the sum of the leftmost `x + 1` elements.
 * @param x The number of elements to consider.
 * @return The calculated difference.
 */
LongType CalculateDifference(int x) {
    return (QueryRight(1, x) - QueryLeft(1, x + 1));
}

/**
 * Processes the query and computes the result based on the current state of the segment tree.
 * @return The result of the query.
 */
LongType ProcessQuery() {
    if (array_size <= 1)
        return 0;
    int left_pointer = 1, right_pointer = array_size - 1, mid_point = array_size / 2;
    if (CalculateDifference(mid_point) <= 0)
        return segment_tree[1].right_sum - segment_tree[1].left_sum;

    // Binary search to find the optimal left and right pointers
    int left = 1, right = mid_point, mid, position;
    if (CalculateDifference(left_pointer) < 0) {
        while (left <= right) {
            mid = (left + right) >> 1;
            if (CalculateDifference(mid) >= 0) {
                position = mid;
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        left_pointer = position;
    }
    if (CalculateDifference(right_pointer) < 0) {
        left = mid_point;
        right = array_size - 1;
        while (left <= right) {
            mid = (left + right) >> 1;
            if (CalculateDifference(mid) >= 0) {
                position = mid;
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        right_pointer = position;
    }

    // Compute the final result
    LongType result = QueryRight(1, right_pointer + 1) - QueryLeft(1, left_pointer);
    if (left_pointer > 1)
        result += QuerySumRight(1, left_pointer - 1) - QuerySumLeft(1, left_pointer - 1);
    if (right_pointer < array_size - 1)
        result += (segment_tree[1].right_sum - QuerySumRight(1, right_pointer + 1)) - (segment_tree[1].left_sum - QuerySumLeft(1, right_pointer + 1));
    return result;
}

/**
 * Solves the problem by processing input and handling queries.
 */
void Solve() {
    cin >> array_size >> query_count;
    for (int i = 1; i <= array_size; i++) {
        cin >> array_values[i];
        unique_values.push_back(array_values[i]);
    }
    for (int i = 1; i <= query_count; i++) {
        int operation_type;
        LongType value;
        cin >> operation_type >> value;
        unique_values.push_back(value);
        queries[i] = {operation_type, 0, value};
    }

    // Sort and deduplicate unique values
    sort(unique_values.begin(), unique_values.end());
    unique_values.erase(unique(unique_values.begin(), unique_values.end()), unique_values.end());
    unique_value_count = unique_values.size();

    // Map values to their IDs and initialize the segment tree
    for (int i = 1; i <= array_size; i++)
        mapped_ids[i] = GetMappedId(array_values[i]);
    for (int i = 1; i <= query_count; i++)
        queries[i].id = GetMappedId(queries[i].value);
    for (int i = 1; i <= array_size; i++)
        UpdateSegmentTree(1, 1, unique_value_count, mapped_ids[i], array_values[i]);

    // Process the initial query
    cout << ProcessQuery() << "\n";

    // Handle each query
    for (int i = 1; i <= query_count; i++) {
        if (queries[i].operation_type == 1) {
            array_size++;
            UpdateSegmentTree(1, 1, unique_value_count, queries[i].id, queries[i].value);
        } else {
            array_size--;
            UpdateSegmentTree(1, 1, unique_value_count, queries[i].id, -queries[i].value);
        }
        cout << ProcessQuery() << "\n";
    }
}

int main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    int test_cases = 1;
    // cin >> test_cases;
    while (test_cases--) {
        Solve();
    }
    return 0;
}
