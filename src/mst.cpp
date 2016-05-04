// Contains the implementations of four different minimum spanning tree 
// algorithms

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits.h>
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
    bool operator<(const Edge& rhs) const {
        return weight < rhs.weight;
    }
};

typedef std::vector<Edge> EdgeContainer;
typedef std::vector<int> VectorSet;

// Data structures for timekeeping
typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds ns;

void ReadInput(char* filename, AdjacencyList &list, AdjacencyMatrix &matrix, int &number_of_vertices, int &number_of_edges);
void SaveResults(char* filename, int v, int e, ns plist, ns pmatrix, ns klist, ns kmatrix);
ns Prim(AdjacencyList &list);
ns Prim(AdjacencyMatrix &matrix);
ns Kruskal(AdjacencyList &list);
ns Kruskal(AdjacencyMatrix &matrix);

void adjust_list(AdjacencyList &list, int v);
void adjust_matrix(AdjacencyMatrix &matrix, int v);
void add_to_list(AdjacencyList &list, int v1, int v2, int weight);
void add_to_matrix(AdjacencyMatrix &matrix, int v1, int v2, int weight);

void read_edges_from_list(AdjacencyList &list, EdgeContainer &container);
void read_edges_from_matrix(AdjacencyMatrix &matrix, EdgeContainer &container);
void sort_edge_container(EdgeContainer &container);
void initialize_set(VectorSet &set, unsigned int v);
bool are_sets_disjoint(VectorSet &set, int v1, int v2);
void join_sets(VectorSet &set, int v1, int v2);

void print_list(AdjacencyList &list);
void print_matrix(AdjacencyMatrix &matrix);
void print_edge_container(EdgeContainer &container);
void print_vector_set(VectorSet &set);


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Please provide 2 filenames." << std::endl;
        return -1;
    }
    
    AdjacencyList list;
    AdjacencyMatrix matrix;
    int number_of_vertices, number_of_edges;
    
    ReadInput(argv[1], list, matrix, number_of_vertices, number_of_edges);

    ns t_prim_list = Prim(list);
    ns t_prim_matrix = Prim(matrix);
    ns t_kruskal_list = Kruskal(list);
    ns t_kruskal_matrix = Kruskal(matrix);


    SaveResults(argv[2], number_of_vertices, number_of_edges, 
            t_prim_list, t_prim_matrix, 
            t_kruskal_list, t_kruskal_matrix);

    
    print_matrix(matrix);
    print_list(list);
}


// Input handler
void ReadInput(char* filename, AdjacencyList &list, AdjacencyMatrix &matrix, int &number_of_vertices, int &number_of_edges) {
    std::ifstream input_file;
    input_file.open(filename);
    
    // Read in parameters
    input_file >> number_of_vertices >> number_of_edges;
    
    // Scale adj_list and adj_matrix
    adjust_list(list, number_of_vertices);
    adjust_matrix(matrix, number_of_vertices);
    
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
void SaveResults(char* filename, int v, int e, ns plist, ns pmatrix, ns klist, ns kmatrix) {
    std::ofstream output_file;
    output_file.open(filename, std::ios::out | std::ios::app );

    output_file << v << ", " << e << ", " 
        << plist.count() << ", " << pmatrix.count() << ", " 
        << klist.count() << ", " << kmatrix.count() << std::endl;

    output_file.close();
}


// Prim with an adjacency list 
ns Prim(AdjacencyList &list) {
    // start clock
    Clock::time_point start = Clock::now();

    EdgeContainer mst;
    std::vector<Node> keys;
    Node tmp;
    tmp.weight = -1;
    tmp.destination = -1;
    keys.push_back(tmp);
    for (unsigned int i = 1; i < list.size(); i++) {
        Node tmp;
        tmp.weight = INT_MAX;
        tmp.destination = -1;
        keys.push_back(tmp);
    }
    
    for(unsigned int i = 1; i < list.size(); i++) {
        // Loop through keys looking for vertices in the mst (==-1)
        for(unsigned int j = 0; j < keys.size(); j++) {
            // node j is already in the mst
            if(keys[j].weight == -1) { 
                for (auto k = list[j].begin(); k != list[j].end(); k++) {
                    int dist_to_mst = k -> weight;
                    if(dist_to_mst > 0 && dist_to_mst < keys[k->destination].weight) {
                        keys[k->destination].weight = dist_to_mst;
                        keys[k->destination].destination = j;
                    }
                }
            }
        }
        
        // Determine the vertex closest to the mst
        int closest_dist = INT_MAX;
        int closest_v = -1;
        for(unsigned int j = 0; j < keys.size(); j++) {
            if(keys[j].weight > 0 && keys[j].weight < closest_dist) {
                closest_dist = keys[j].weight;
                closest_v = j;
            }
        }
        Edge e;
        e.v1 = closest_v;
        e.v2 = keys[closest_v].destination;
        e.weight = keys[closest_v].weight;
        mst.push_back(e);

        // Set closest vertex to -1 to "include" it in the mst
        keys[closest_v].weight = -1;
    }

    //printf("Minimum Spanning Tree Edges:\n");
    //for(auto e : mst) printf("%d->%d of weight %d\n", e.v1, e.v2, e.weight);
    //printf("\n");

    // stop clock and return time
    Clock::time_point end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}


// Prim with an adjacency matrix
ns Prim(AdjacencyMatrix &matrix) {
    // start clock
    Clock::time_point start = Clock::now();

    EdgeContainer mst;
    std::vector<Node> keys;
    Node tmp;
    tmp.weight = -1;
    tmp.destination = -1;
    keys.push_back(tmp);
    for (unsigned int i = 1; i < matrix.size(); i++) {
        Node tmp;
        tmp.weight = INT_MAX;
        tmp.destination = -1;
        keys.push_back(tmp);
    }
    
    for(unsigned int i = 1; i < matrix.size(); i++) {
        // Loop through keys looking for vertices in the mst (==-1)
        for(unsigned int j = 0; j < keys.size(); j++) {
            // node j is already in the mst
            if(keys[j].weight == -1) { 
                for(unsigned int k = 0; k < matrix[j].size(); k++) {
                    int dist_to_mst = matrix[j][k];
                    if(dist_to_mst > 0 && dist_to_mst < keys[k].weight) {
                        keys[k].weight = matrix[j][k];
                        keys[k].destination = j;
                    }
                }
            }
        }
        
        // Determine the vertex closest to the mst
        int closest_dist = INT_MAX;
        int closest_v = -1;
        for(unsigned int j = 0; j < keys.size(); j++) {
            if(keys[j].weight > 0 && keys[j].weight < closest_dist) {
                closest_dist = keys[j].weight;
                closest_v = j;
            }
        }
        Edge e;
        e.v1 = closest_v;
        e.v2 = keys[closest_v].destination;
        e.weight = keys[closest_v].weight;
        mst.push_back(e);

        // Set closest vertex to -1 to "include" it in the mst
        keys[closest_v].weight = -1;
    }

    //printf("Minimum Spanning Tree Edges:\n");
    //for(auto e : mst) printf("%d->%d of weight %d\n", e.v1, e.v2, e.weight);
    //printf("\n");

    // stop clock and return time
    Clock::time_point end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}


// Kruskal with an adjacency list 
ns Kruskal(AdjacencyList &list) {
    // start clock
    Clock::time_point start = Clock::now();
    
    // Organize all edges, smallest weight first
    EdgeContainer all_edges;
    read_edges_from_list(list, all_edges);
    sort_edge_container(all_edges);

    // Place each vector in it's own disjoint set
    VectorSet set;
    initialize_set(set, list.size());

    EdgeContainer mst; // Will contain the edges of minimum spanning tree

    for (Edge e : all_edges) {
        if (are_sets_disjoint(set, e.v1, e.v2)) {
            mst.push_back(e);
            join_sets(set, set[e.v1], set[e.v2]);
        }
    }

    // stop clock and return time
    Clock::time_point end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}


// Kruskal with an adjacency matrix
ns Kruskal(AdjacencyMatrix &matrix) {
    // start clock
    Clock::time_point start = Clock::now();

    // Organize all edges, smallest weight first
    EdgeContainer all_edges;
    read_edges_from_matrix(matrix, all_edges);
    sort_edge_container(all_edges);

    // Place each vector in it's own disjoint set
    VectorSet set;
    initialize_set(set, matrix.size());

    EdgeContainer mst; // Will contain the edges of minimum spanning tree

    for (Edge e : all_edges) {
        if (are_sets_disjoint(set, e.v1, e.v2)) {
            mst.push_back(e);
            join_sets(set, set[e.v1], set[e.v2]);
        }
    }

    // stop clock and return time
    Clock::time_point end = Clock::now();
    return std::chrono::duration_cast<ns> (end - start);
}



// Helper Functions
void adjust_list(AdjacencyList &list, int v) {
    list.resize(v);
}

void adjust_matrix(AdjacencyMatrix &matrix, int v) {
    matrix.resize(v);

    for (int i = 0; i < v; i++) {
        matrix[i].resize(v,-1);
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

void read_edges_from_list(AdjacencyList &list, EdgeContainer &container) {
    /* This code does not check if an edge has already been included.
     * Since these are undirected graphs each edge will be included twice.
     * Kruskal's algorithm checks if the vectors of the edge are disjointed, 
     * so repeat edges will not affect the results.
     * It would take n^2 time to check if an edge has already been included,
     * meanwhile the cost of redundant edges is 2n. I decided to use the
     * solution with the least bottleneck, since we are coding for time.
     * */
    for (unsigned int i = 0; i < list.size(); i++) {
        for (auto n = list[i].begin(); n != list[i].end(); n++) {
            int v1 = i;
            int v2 = n -> destination;
            if (v1 > v2)
                std::swap(v1,v2);

            Edge e;
            e.v1 = v1;
            e.v2 = v2;
            e.weight = n -> weight;
            container.push_back(e);
        }
    }
}

void read_edges_from_matrix(AdjacencyMatrix &matrix, EdgeContainer &container) {
    for (unsigned int i = 0; i < matrix.size(); i++) {
        for (unsigned int j = 0; j < i; j++) {
            if (matrix[i][j] == -1)
                continue;

            Edge e;
            e.v1 = j;
            e.v2 = i;
            e.weight = matrix[i][j];
            container.push_back(e);
        }
    }
}

void sort_edge_container(EdgeContainer &container) {
    std::sort(container.begin(), container.end());
}
void initialize_set(VectorSet &set, unsigned int v) {
    set.resize(v); // each index represents a vector
    for (unsigned int i = 0; i < v; i++) {
        set[i] = i; // set each vector to its own index, represents unique sets
    }
}

bool are_sets_disjoint(VectorSet &set, int v1, int v2) {
    return ((set[v1] != set[v2])? true : false);
}

void join_sets(VectorSet &set, int a, int b) {
    int conquering_set = (a < b) ? a : b;
    int merging_set = (a < b) ? a: b;

    for (unsigned int i = 0; i < set.size(); i++) {
        if (set[i] == merging_set)
            set[i] = conquering_set;
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
            if (matrix[i][j] == -1)
                std::cout << "*\t";
            else
                std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_edge_container(EdgeContainer &container) {
    std::cout << "Edge Container" << std::endl;

    for (unsigned int i = 0; i < container.size(); i++) {
        std::cout << i << ".\t" <<  container[i].v1 << "->" << container[i].v2
            << ", w: " << container[i].weight << std::endl;

    }
    std::cout << std::endl;
}

void print_vector_set(VectorSet &set) {
    std::cout << "Vector Set" << std::endl;

    for (unsigned int i = 0; i < set.size(); i++) {
        std::cout << i << " [" << set[i] << "] -- ";
    }
    std::cout << std::endl << std::endl;
}
