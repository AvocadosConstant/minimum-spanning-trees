// Contains the implementations of four different minimum spanning tree 
// algorithms

#include <iostream>
#include <fstream>
#include <vector>
#include <forward_list>

struct Node {
    int destination;
    int weight;
};

typedef std::vector< std::forward_list<Node> > AdjacencyList;
typedef std::vector< std::vector<int> > AdjacencyMatrix;


void ReadInput(char* filename, AdjacencyList &list, AdjacencyMatrix &matrix);

void adjust_list(AdjacencyList &list, int v);
void adjust_matrix(AdjacencyMatrix &matrix, int v);

void add_to_list(AdjacencyList &list, int v1, int v2, int weight);
void add_to_matrix(AdjacencyMatrix &matrix, int v1, int v2, int weight);

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

    print_list(list);
    print_matrix(matrix);
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
// Prim with an adjacency list 
// Prim with an adjacency matrix
// Kruskal with an adjacency list 
// Kruskal with an adjacency matrix


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

void print_list(AdjacencyList &list) {
    std::cout << "Adjacency List" << std::endl;

    for (int i = 0; i < list.size(); i++) {
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
    for (int i = 0; i < matrix.size(); i++) {
        std::cout << i << "\t";
    }
    std::cout << std::endl;

    for (int i = 0; i < matrix.size(); i++) {
        std::cout << i << "\t";
        for (int j = 0; j < matrix[i].size(); j++) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

