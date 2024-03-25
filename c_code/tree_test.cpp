#include "tree_test.h"
#include <cassert>

using namespace std;
std::vector<bool> results;

void customAssert(string expected, string actual, const std::string& message) {
    if (!(expected==actual)) {
        std::cout << "Assertion failed: " << message << std::endl;
        cout << "expected: " << expected << endl;
        cout << "actual: " << actual << endl;
        results.push_back(false);
    } else {
        std::cout << "Assertion passed: " << message << std::endl;
        results.push_back(true);
    }
}



int main() {
    cout << "testing printNewick and pruneExtinctSpecies" << endl;

    // test case 1
    speciation_tree_node* root1 = createNode(1, 0);
    root1-> left_child = createNode(1,10);
    root1-> right_child = createNode(2,10);
    root1->left_child->left_child = createNode(1,20);
    root1->left_child->right_child = createNode(3,20);
    root1->right_child->left_child = createNode(2,30);
    root1->right_child->right_child = createNode(4,30);
    set<cell_type> extant_species1({2,3,4});
    string expected1 = "((1:10,3:10):10,(2:20,4:20):10):0";
    string actual1 = printNewick(root1);
    string expected1p = "((3:10):10,(2:20,4:20):10):0";
    string actual1p = printNewick(pruneExtinctSpecies(root1,extant_species1));
    customAssert(expected1,actual1, "test case 1 print");
    customAssert(expected1p,actual1p, "test case 1 prune");

    // test case 2
    speciation_tree_node* root2 = createNode(1, 0);
    root2-> left_child = createNode(1,10);
    root2-> right_child = createNode(2,10);
    root2->left_child->left_child = createNode(1,20);
    root2->left_child->right_child = createNode(3,20);
    root2->right_child->left_child = createNode(2,30);
    root2->right_child->right_child = createNode(4,30);
    set<cell_type> extant_species2({1,2,4});
    string expected2 = "((1:10,3:10):10,(2:20,4:20):10):0";
    string actual2 = printNewick(root2);
    string expected2p = "(1:20,(2:20,4:20):10):0";
    string actual2p = printNewick(pruneExtinctSpecies(root2,extant_species2));
    customAssert(expected2,actual2, "test case 2 print");
    customAssert(expected2p,actual2p, "test case 2 prune");

    // test case 3
    speciation_tree_node* root3 = createNode(1, 0);
    root3-> left_child = createNode(1,10);
    root3-> right_child = createNode(2,10);

    root3->left_child->left_child = createNode(1,20);
    root3->left_child->right_child = createNode(3,20);
    root3->right_child->left_child = createNode(2,30);
    root3->right_child->right_child = createNode(4,30);

    root3->left_child->left_child->left_child = createNode(1,40);
    root3->left_child->left_child->right_child = createNode(5,40);
    root3->left_child->right_child->left_child = createNode(3,50);
    root3->left_child->right_child->right_child = createNode(6,50);
    set<cell_type> extant_species3({2,4});
    string expected3 = "(((1:20,5:20):10,(3:30,6:30):10):10,(2:20,4:20):10):0";
    string actual3 = printNewick(root3);
    string expected3p = "((2:20,4:20):10):0";
    string actual3p = printNewick(pruneExtinctSpecies(root3,extant_species3));
    customAssert(expected3,actual3, "test case 3 print");
    customAssert(expected3p,actual3p, "test case 3 prune");

    for (const auto& result : results) {
        std::cout << (result ? "\033[1;32m+\033[0m" : "\033[1;31m*\033[0m");
    }
    std::cout << endl;
    return 0;
}