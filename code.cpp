#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <chrono>

// The vertex structure
struct vertex
{
    int w;
    std::vector<int> adjEdges;
};

// The vector that holds the vertices
std::vector<vertex> vertices;

// The vector that holds a (ordered) pair of vertices for each edge
std::vector<std::pair<int,int>> edges;

// Our dynamic programming variables (for the edges)
std::vector<int> sumWeights; // [i] = W(i, edges[i].second)
std::vector<int> firstResult; // [i] = S'(i, edges[i].second)
std::vector<int> biggestWeightEdge; // [i] = argmax W

std::vector<int> bestResult; // [i] = min[x] S(i, edges[i].second)(x)
std::vector<int> bestResultEdge; // [i] = the edge which second is the best result of i
std::vector<int> bestResultDist; // [i] = the distance from i to the best result of i
// basically the "inverse" of biggestWeightEdge: fromPath[i] = e <=> biggestWeightEdge[e] = i
std::vector<int> fromPath;

// Utility functions
void computeSupportValuesForEdge(int e);
void solveMinForEdge(int e);

int main()
{
    // Read the number of vertices
    int n;
    std::cin >> n;

    std::cout << "n = " << n << std::endl;
    edges.resize(2*(n-1)); // Each edge is stored with its opposite
    vertices.resize(n); // n vertices

    // For each of the n-1 edges, append them to the 
    for (int i = 0; i < n-1; i++)
    {
        int u, v;
        std::cin >> u >> v;

        // Append the edges to the list, since we need them
        edges[2*i] = std::make_pair(u-1, v-1);
        edges[2*i+1] = std::make_pair(v-1, u-1);

        // Append the adjacent edges to each vertex list
        vertices[u-1].adjEdges.push_back(2*i);
        vertices[v-1].adjEdges.push_back(2*i+1);
    }

    // For each of the vertices, find its weight
    for (int i = 0; i < n; i++)
        std::cin >> vertices[i].w;

    // Resize the support functions
    sumWeights.resize(2*(n-1), -1);
    firstResult.resize(2*(n-1), -1);
    biggestWeightEdge.resize(2*(n-1), -1);
    bestResult.resize(2*(n-1), -1);
    bestResultEdge.resize(2*(n-1), -1);
    bestResultDist.resize(2*(n-1), -1);
    fromPath.resize(2*(n-1), -1);

    // The variable to hold the min for each edge
    int smin = std::numeric_limits<int>::max();

    // Time from here
    auto then = std::chrono::high_resolution_clock::now();

    // For each of the edges of the list, we compute its support functions
    for (int i = 0; i < 2*(n-1); i++)
        computeSupportValuesForEdge(i);

    // Now, we compute the values and edges for their twins
    for (int i = 0; i < n-1; i++)
    {
        solveMinForEdge(2*i);
        solveMinForEdge(2*i+1);
        smin = std::min(smin, bestResult[2*i] + bestResult[2*i+1]);
    }

    // Time until here
    auto now = std::chrono::high_resolution_clock::now();

    // Output the result
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(now - then).count();
    std::cout << "result = " << smin << std::endl;
    std::cout << "time = " << us << "us" << std::endl;
}

// Here, we compute the values of W and S' for each edge
void computeSupportValuesForEdge(int e)
{
    // Bail out if the result was already computed
    if (sumWeights[e] != -1) return;

    // The "implicit" coding is that the edge is i
    // and u here will be edges[e].second (why? because, for
    // all edges coded in vertices[v].adjEdges, u = edges[...].first!)
    int u = edges[e].second;

    sumWeights[e] = vertices[u].w;
    firstResult[e] = 0;

    for (int en : vertices[u].adjEdges)
    {
        // Ignore the same edge, important
        if (edges[en].second == edges[e].first) continue;

        // Recursively compute the support values and add them
        computeSupportValuesForEdge(en);
        sumWeights[e] += sumWeights[en];
        firstResult[e] += sumWeights[en] + firstResult[en];

        if (biggestWeightEdge[e] == -1
            || sumWeights[biggestWeightEdge[e]] < sumWeights[en])
            biggestWeightEdge[e] = en;
    }

    // Compute fromPath from here
    if (biggestWeightEdge[e] != -1)
        fromPath[biggestWeightEdge[e]] = e;
}

// Here, we will actually do the tree traversal to find the minimum s
void solveMinForEdge(int e)
{
    // If the result has already been solved, return
    if (bestResult[e] != -1) return;

    // Assume the first result is the best result for now
    bestResult[e] = firstResult[e];
    bestResultEdge[e] = e;
    bestResultDist[e] = 0;

    // If there are no subtrees in that tree, the best result is what we had already found
    int biggest = biggestWeightEdge[e];
    if (biggest == -1) return;

    // Compute the mininmum value for the biggest edge
    solveMinForEdge(biggest);

    // And travel backwards for the solution
    int curDist = 1 + bestResultDist[biggest];
    int curEdge = bestResultEdge[biggest];
    
    // The current result is computed using the formula described above
    int curResult = curDist * (sumWeights[e] - sumWeights[biggest])
        + firstResult[e] - firstResult[biggest] - sumWeights[biggest]
        + bestResult[biggest];

    // Travel backwards searching for all solutions
    while (curDist > 0)
    {
        // Change the result
        bestResult[e] = curResult;
        bestResultEdge[e] = curEdge;
        bestResultDist[e] = curDist;

        // And stop as soon as the result subtraction turns negative
        if (sumWeights[e] - 2 * sumWeights[curEdge] < 0) break;
        curResult -= sumWeights[e] - 2 * sumWeights[curEdge];
        curEdge = fromPath[curEdge];
        curDist--;
    }
}