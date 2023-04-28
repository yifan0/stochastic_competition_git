#include "tree.cpp"
#include <iostream> // std::cout
using namespace std;

int main() {
    // start with a tree the same as the simulation
    set<speciation_tree_node*> species;
    speciation_tree_node* speciation_root = new speciation_tree_node(1, 0, nullptr); // root has value 1, the starting value
    speciation_tree_node* starting_value = new speciation_tree_node(0.7, 0, speciation_root);
    speciation_root->left_child = starting_value;
    species.insert(starting_value);

    speciation_tree_node* A = new speciation_tree_node(0.5, 1, nullptr);
    speciation_event(starting_value, A);

    cout << toString(speciation_root) << endl;
    
    speciation_tree_node* B = new speciation_tree_node(0.1, 2, nullptr);
    speciation_event(A, B);
    species.insert(B);

    cout << toString(speciation_root) << endl;

    speciation_tree_node* C = new speciation_tree_node(0.2, 3, nullptr);
    speciation_event(A, C);
    species.insert(C);

    cout << toString(speciation_root) << endl;

    speciation_tree_node* D = new speciation_tree_node(0.3, 4, nullptr);
    speciation_event(C, D);
    species.insert(D);

    speciation_tree_node* E = new speciation_tree_node(0.4, 5, nullptr);
    speciation_event(B, E);
    species.insert(E);

    cout << toString(speciation_root) << endl;
    
    println("Depth of tree before pruning = %d", get_depth(speciation_root)); 
    prune(speciation_root, species);
    println("Depth of tree after pruning = %d", get_depth(speciation_root)); 
    cout << toString(speciation_root) << endl;
    
    return 0;
}
