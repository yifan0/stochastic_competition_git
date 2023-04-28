#include <cstdio>
#include <string.h>
#include <tuple>
#include <set>
using namespace std;

#ifndef TREE_

#define println(...) { printf(__VA_ARGS__); printf("\n"); }

typedef double cell_type;

struct speciation_tree_node {
    cell_type val;
    size_t time;
    speciation_tree_node* parent;
    speciation_tree_node* right_child;
    speciation_tree_node* left_child;
    speciation_tree_node(cell_type v, size_t t, speciation_tree_node* p) : val(v), time(t), parent(p) {
        right_child = nullptr;
        left_child = nullptr;
    }
};

typedef std::tuple<int, int, speciation_tree_node*> cell_update;

speciation_tree_node* find_node(speciation_tree_node* root, cell_type target) {
    if(root->val == target) {
        return root;
    }

    speciation_tree_node* result = find_node(root->right_child, target);
    if(result != nullptr) return result;
    
    result = find_node(root->left_child, target);
    if(result != nullptr) return result;

    return nullptr;
}

// caller is responsible for adding a semicolon at the end
string toString(speciation_tree_node* root) {
    // for a leaf, just return the label
    if(root->right_child == nullptr || root->left_child == nullptr) {
        return std::to_string(root->val);
    }

    // for a parent, group the children in ()
    string output = "(" + toString(root->right_child) + "," + toString(root->left_child) + ")";
    //output += std::to_string(root->val); // do not label internal nodes
    return output;
}

void delete_tree(speciation_tree_node* root) {
    if(root == nullptr) return;
    delete_tree(root->left_child);
    delete_tree(root->right_child);
    delete root;
}

size_t get_depth(speciation_tree_node* root) {
    if(root == nullptr) return 0;
    size_t depth = 0;
    depth = max(depth, get_depth(root->left_child));
    depth = max(depth, get_depth(root->right_child));
    return depth+1;
}

void speciation_event(speciation_tree_node* old_species, speciation_tree_node* new_species) {

    if(old_species->parent->left_child == nullptr) {
        old_species->parent->left_child = new_species;
        new_species->parent = old_species->parent;
        return;
    } else if(old_species->parent->right_child == nullptr) {
        old_species->parent->right_child = new_species;
        new_species->parent = old_species->parent;
        return;
    }

    // create new node for parent
    speciation_tree_node* internal_node = new speciation_tree_node(old_species->val, old_species->time, old_species->parent);
    internal_node->left_child = old_species;
    internal_node->right_child = new_species;

    // insert new node for new species
    if(old_species->parent->left_child == old_species) {
        old_species->parent->left_child = internal_node;
    } else if(old_species->parent->right_child == old_species) {
        old_species->parent->right_child = internal_node;
    } else {
        println("error adding new node");
    }

    new_species->parent = internal_node;
    old_species->parent = internal_node;
}

void normalize(speciation_tree_node* root, cell_type avg) {
    if(root == nullptr) return;
    root->val = root->val/avg;
    normalize(root->left_child, avg);
    normalize(root->right_child, avg);
}

void prune(speciation_tree_node* root, set<speciation_tree_node*> species) {
    if(root == nullptr) return;
    if(species.count(root) == 0 && root->left_child == nullptr && root->right_child == nullptr) {
        // species is extinct and so should be removed
        speciation_tree_node* parent = root->parent;
        speciation_tree_node* sibling;
        if(root->parent->left_child == root) {
            sibling = parent->right_child;
            parent->right_child = nullptr;
        }
        else if(root->parent->right_child == root) {
            sibling = parent->left_child;
            parent->left_child = nullptr;
        }
        else {
            println("Error pruning species");
        }

        if(parent->parent->left_child == parent) {
            parent->parent->left_child = sibling;
        } else {
            parent->parent->right_child = sibling;
        }
        delete(root);
        delete(parent);

    }
    prune(root->left_child, species);
    prune(root->right_child, species);
}

#endif
