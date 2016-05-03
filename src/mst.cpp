// Contains the implementations of four different minimum spanning tree 
// algorithms

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <forward_list>
#include <chrono>

// Data structures for graph representations
struct Node {
    int destination;
    int weight;
};

typedef std::vector< std::forward_list<Node> > AdjacencyList;
typedef std::vector< std::vector<int> > AdjacencyMatrix;

// Data structures for Kruskal's algorithms
struct Edge {
    int v1;
    int v2;
    int weight;
};

typedef std::vector<Edge> EdgeContainer;
typedef std::vector<int> VectorSet;

// Data structures for timekeeping
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds ns;

void ReadInput(char* filename, AdjacencyList &list, AdjacencyMatrix &matrix);
void SaveResults(char* filename, ns plist, ns pmatrix, ns klist, ns kmatrix);
ns Prim(AdjacencyList &list);
ns Prim(AdjacencyMatrix &matrix);
ns Kruskal(AdjacencyList &list);
ns Kruskal(AdjacencyMatrix &matrix);

void adjust_list(AdjacencyList &list, int v);
void adjust_matrix(AdjacencyMatrix &matrix, int v);
void add_to_list(AdjacencyList &list, int v1, int v2, int weight);
void add_to_matrix(AdjacencyMatrix &matrix, int v1, int v2, int weight);
void read_edges_from_list(AdjacencyList &list, EdgeContainer &all_edges);
void print_list(AdjacencyList &list);
void print_matrix(AdjacencyMatrix &matrix);


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Please provide 2 filenames." << std::endl;
        return -1;
    }
    
    AdjacencyList list;
    AdjacencyMatrix matrix;
    
    ReadInput(argv[1], list, matrix);

    ns t_prim_list = Prim(list);
    ns t_prim_matrix = Prim(matrix);
    ns t_kruskal_list = Kruskal(list);
    ns t_kruskal_matrix = Kruskal(matrix);

    SaveResults(argv[2], t_prim_list, t_prim_matrix, 
            t_kruskal_list, t_kruskal_matrix);

    
    print_matrix(matrix);
    print_list(list);
}


// Input handler
void ReadInput(char* filename, AdjacencyList &list, AdjacencyMatrix &matrix) {
    std::ifstream input_file;
    input_file.open(filename);
    
    // Read in parameters
    int number_of_vectors, number_of_edges;
    input_file >> number_of_vectors >> number_of_edges;
    
    // Scale adj_list and adj_matrix
    adjust_list(list, number_of_vectors);
    adjust_matrix(matrix, number_of_vectors);
    
    // Read in Edges
    int v1, v2, weight;
    for (int i = 0; i < number_of_edges; i++) {
        input_file >> v1 >> v2 >> weight;

        add_to_list(list, v1, v2, weight);
        add_to_matrix(matrix, v1, v2, weight);
    }
    input_file.close();
}

// Output handler
void SaveResults(char* filename, ns plist, ns pmatrix, ns klist, ns kmatrix) {
    std::ofstream output_file;
    output_file.open(filename);

    output_file << plist.count() << ", " << pmatrix.count() << ", " 
        << klist.count() << ", " << kmatrix.count() << std::endl;

    output_file.close();
}

// Prim with an adjacency list 
ns Prim(AdjacencyList &list) {
    // start clock
    Clock::time_point start, end;
    start = Clock::now();

    /* PSEUDOCODE
     *
     * Q = V
     * key[v] = infinity for all v in V
     * key[s] = 0 // s is the starting vertex
     *
     * while (!Q.empty()):
     *      u = Q.extract_min()
     *      for each v adjacent to u:
     *          if (v in Q and u->v.weight < key[v]):
     *              key[v] = u->v.weight
     *              pi[v] = u
     */

    // stop clock and return time
    end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}

// Prim with an adjacency matrix
ns Prim(AdjacencyMatrix &matrix) {
    // start clock
    Clock::time_point start, end;
    start = Clock::now();

    /* PSEUDOCODE
     *
     * Q = V
     * key[v] = infinity for all v in V
     * key[s] = 0 // s is the starting vertex
     *
     * while (!Q.empty()):
     *      u = Q.extract_min()
     *      for each v adjacent to u:
     *          if (v in Q and u->v.weight < key[v]):
     *              key[v] = u->v.weight
     *              pi[v] = u
     */

    // stop clock and return time
    end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}

// Kruskal with an adjacency list 
ns Kruskal(AdjacencyList &list) {
    // start clock
    Clock::time_point start, end;
    start = Clock::now();

    /* PSEUDOCODE
     *
     * A = null // A is a list of edges. stores MST
     *
     * form a set out of each vertex
     * sort the edges into nondecreasing order by weight
     * loop through each e in sorted edges:
     *      if the two sets u,v (vertices of e) are disjoint:
     *          add the edge to A
     *          unite the two sets u + v
    */
    
    EdgeContainer all_edges;
    read_edges_from_list(list, all_edges);
    // Sort Edges, smallest first
    // Print EdgeContainer

    VectorSet set;
    // Resize set to number of vectors
    // Initialize set

    EdgeContainer mst;

    for (Edge e : all_edges) {
        // compare v1 and v2 in disjoint_set
        // add edge to mst
        // unite v1 and v2
    }

    // stop clock and return time
    end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}

// Kruskal with an adjacency matrix
ns Kruskal(AdjacencyMatrix &matrix) {
    // start clock
    Clock::time_point start, end;
    start = Clock::now();

    /* PSEUDOCODE
     *
     * A = null // A is a list of edges. stores MST
     *
     * form a set out of each vertex
     * sort the edges into nondecreasing order by weight
     * loop through each e in sorted edges:
     *      if the two sets u,v (vertices of e) are disjoint:
     *          add the edge to A
     *          unite the two sets u + v
    */

    // stop clock and return time
    end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}


// Helper Functions
void adjust_list(AdjacencyList &list, int v) {
    list.resize(v);
}

void adjust_matrix(AdjacencyMatrix &matrix, int v) {
    matrix.resize(v);

    for (int i = 0; i < v; i++) {
        matrix[i].resize(v,0);
    }
}

void add_to_list(AdjacencyList &list, int v1, int v2, int weight) {
    if (v1 == v2) 
        return;
    
    Node n1;
    n1.weight = weight;
    n1.destination = v2;
    list[v1].push_front(n1);

    Node n2;
    n2.weight = weight;
    n2.destination = v1;
    list[v2].push_front(n2);
}

void add_to_matrix(AdjacencyMatrix &matrix, int v1, int v2, int weight) {
    matrix[v1][v2] = weight;
    matrix[v2][v1] = weight;
}

void read_edges_from_list(AdjacencyList &list, EdgeContainer &all_edges) {
    for (int i = 0; i < list.size(); i++) {
        for (auto n = list[i].begin(); n != list[i].end(); n++) {
            int v1 = i;
            int v2 = n -> destination;
            if (v1 > v2)
                std::swap(v1,v2);

            Edge e;
            e.v1 = v1;
            e.v2 = v2;
            e.weight = n -> weight;
            all_edges.push_back(e);
        }
    }
}

void print_list(AdjacencyList &list) {
    std::cout << "Adjacency List" << std::endl;

    for (unsigned int i = 0; i < list.size(); i++) {
        std::cout << i << ":\t";

        for (auto n = list[i].begin(); n != list[i].end(); n++) {
            std::cout << "[" << n->destination << "," << n->weight << "]\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_matrix(AdjacencyMatrix &matrix) {
    std::cout << "Adjacency Matrix" << std::endl;

    std::cout << "\t";
    for (unsigned int i = 0; i < matrix.size(); i++) {
        std::cout << i << "\t";
    }
    std::cout << std::endl;

    for (unsigned int i = 0; i < matrix.size(); i++) {
        std::cout << i << "\t";
        for (unsigned int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

