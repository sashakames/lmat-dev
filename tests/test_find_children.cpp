#include "../include/all_headers.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cassert>

using namespace std;
using namespace metag;

char *base = "/usr/mic/post1/metagenomics/small_test_data/";

int main(int argc, char *argv[]) {
  char buf[1024];

  sprintf(buf, "%s/tax_nodes.dat", base);
  TaxTree tax_tree(buf);
  tax_tree.findChildren();

  for (map<int, TaxNode*>::const_iterator t = tax_tree.m_tree.begin(); t != tax_tree.m_tree.end(); t++) {
    cout << endl;
    t->second->printChildren();
  }
}
