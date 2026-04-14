#include <iostream>
#include <queue>
#include <set>
using namespace std;

typedef int cell_type;

// speciation_tree_node structure for the binary tree
struct speciation_tree_node {
    cell_type val;
    int time;
    speciation_tree_node* left_child;
    speciation_tree_node* right_child;
    speciation_tree_node(cell_type value, int t) : val(value), time(t), left_child(nullptr), right_child(nullptr) {}
};

// Function to create a new speciation_tree_node
speciation_tree_node* createNode(cell_type value, int t) {
    return new speciation_tree_node(value, t);
}




speciation_tree_node* pruneExtinctSpecies(speciation_tree_node* root, const std::set<cell_type>& extant_species) {
    if (root == nullptr) {
        return nullptr;
    }

    root->left_child = pruneExtinctSpecies(root->left_child, extant_species);
    root->right_child = pruneExtinctSpecies(root->right_child, extant_species);
    // prune extinct leaf
    bool isExtinct = extant_species.find(root->val) == extant_species.end();
    if (isExtinct && root->left_child == nullptr && root->right_child == nullptr) {
        delete root;
        return nullptr;
    }
    // prune speciation_tree_nodes with a single child
    if (root->left_child == nullptr && root->right_child != nullptr) {
        // (not prune if the only child is the right child, preserving the mutation history)
        speciation_tree_node* right_childChild = root->right_child;
        return root;
    } else if (root->left_child != nullptr && root->right_child == nullptr) {
        // prune if the only chid is the left child, which has the same value as the root since no mutation occurred.
        speciation_tree_node* left_childChild = root->left_child;
        delete root;
        return left_childChild;
    }

    return root;
}

// Function to print the binary tree in Newick format with distances to parent nodes
string printNewick(speciation_tree_node* root, size_t parentTime = 0) {
    if (root == nullptr) {
        return "";
    }

    if (root->left_child == nullptr && root->right_child == nullptr) {
        return to_string(root->val) + ":" + to_string(root->time - parentTime);
    } else {
        // add comma if both children are not nullptr. If any of the child is nullptr, expect only one child in the parantheses with no comma.
        string addComma = (root->left_child != nullptr && root->right_child != nullptr) ? ",": "";
        return "(" + printNewick(root->left_child, root->time) + addComma + printNewick(root->right_child, root->time) + ")" + ":" + to_string(root->time - parentTime);
    }
}