#include <cstdio>
#include <string>
#include <tuple>
#include <set>
#include <iostream>
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
    // add branch length
    int branch_len_left = root->left_child->time - root->time;
    int branch_len_right = root->right_child->time - root->time;
    string output = "(" + toString(root->right_child) + ":" + to_string(branch_len_left) + "," + toString(root->left_child) + ":" + to_string(branch_len_right) + ")";
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
    // only occurs at the first speciation event
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
    old_species->time = new_species->time;
}

void normalize(speciation_tree_node* root, cell_type avg) {
    if(root == nullptr) return;
    root->val = root->val/avg;
    normalize(root->left_child, avg);
    normalize(root->right_child, avg);
}

void prune(speciation_tree_node* root, set<speciation_tree_node*> species) {
    if(root == nullptr){
        return;
    } 
    if(species.count(root) == 0 && root->left_child == nullptr && root->right_child == nullptr) {
        // return when parent of parent is nullptr
        if (root->parent->parent == nullptr){
            cout << "missed extinct species " << root->val << endl;
            return;
        }
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
            sibling->parent = parent->parent;
        } else {
            parent->parent->right_child = sibling;
            sibling->parent = parent->parent;
        }
        delete(root);
        delete(parent);
        return;
    }
    prune(root->left_child, species);
    prune(root->right_child, species);
}

void prune_test(speciation_tree_node* root, set<speciation_tree_node*> species, speciation_tree_node* speciation_root) {
    if(root == nullptr){
        cout << "root null" << endl;
        return;
    } 
    cout << "root called " << root->val << " at depth " << get_depth(root) << endl;
    cout << "tree at call " << toString(speciation_root).c_str() << endl;
    if(species.count(root) == 0 && root->left_child == nullptr && root->right_child == nullptr) {
        cout << "root to delete " << root->val << endl;
        cout << toString(speciation_root).c_str() << endl;
        if (root->parent->parent == nullptr){
            cout << "parent of parent is null" << endl;
            cout << toString(speciation_root).c_str() << endl;
            return;
        }
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
            sibling->parent = parent->parent;
        } else {
            parent->parent->right_child = sibling;
            sibling->parent = parent->parent;
        }
        cout << "root = " << root->val << " deleted " << endl;
        delete(root);
        delete(parent);
        // cout << "tree after deletion " << toString(speciation_root).c_str() << endl;
        return;
    }
    // left child
    if (root->left_child == nullptr){
        cout << "left child of " << root->val << " at depth " << get_depth(root) << " is null" << endl;
    } else {
        cout << "left child of " << root->val << " at depth " << get_depth(root) << " is " << root->left_child->val << " at depth " << get_depth(root->left_child) << endl;
    }
    prune_test(root->left_child, species, speciation_root);
    // right child
    if (root->right_child == nullptr){
        cout << "right child of " << root->val << " at depth " << get_depth(root) << " is null" << endl;
    } else {
        cout << "right child of " << root->val << " at depth " << get_depth(root) << " is " << root->right_child->val << " at depth " << get_depth(root->right_child) << endl;
    }
    prune_test(root->right_child, species, speciation_root);
}

#endif
