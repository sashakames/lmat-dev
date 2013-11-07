#include <iostream>
#include <string>
#include <omp.h>
#include "all_headers.hpp"

using namespace std;
using namespace metag;


int main(int argc, char *argv[]) {

int nthreads;
#pragma omp parallel
{
  int num = omp_get_num_threads();
  int id = omp_get_thread_num();
  if (id == 0) {
    nthreads = num;
  }
}

cout << "num threads: " << nthreads << endl;

  TaxTree tax_tree("../../../../runtime_inputs/kpath_taxonomy.dat", "../../../../runtime_inputs/gid_to_kpath_tax_id_microbe2.dat", nthreads);

  set<uint32_t> tax_ids;
  tax_ids.insert(974);
  tax_ids.insert(13910);
  tax_ids.insert(48497);
  tax_ids.insert(157975);
  tax_ids.insert(157976);
  tax_ids.insert(157999);
  tax_ids.insert(190173);
  tax_ids.insert(273213);
  tax_ids.insert(295518);
  tax_ids.insert(539837);

  set<uint32_t> lowest;

cout << "running threaded test:\n";
#pragma omp parallel
{
  omp_set_num_threads(5);
  int id = omp_get_thread_num();
  int num = omp_get_num_threads();
  if (id == 0) {
    cout << "num threads: " << num << endl << endl;
  }
  vector<uint32_t> r;
  for (int j=0; j<20; j++) {
    const TaxNode *n = tax_tree.lowestCommon(tax_ids, id);
    if (n == 0) {
      r.push_back(1000000000);
    } else {
      r.push_back(n->id());
    }
  }
  #pragma omp critical
  {
    cout << "id: " << id << " results: ";
    for (int h=0; h<20; h++) cout << r[h] << " ";
    cout << endl;
  }
}

  cout << "calling TaxTree::lowestCommon(set<uint32_t>)\n";
  const TaxNode *n = tax_tree.lowestCommon(tax_ids);
  cout << "lowestCommon returned: " << n << endl;
  if (n != 0) {
    TaxNode *nn = const_cast<TaxNode*>(n);
    nn->write();
  }

  cout << "\ncalling TaxTree::lowestCommon(set<uint32_t>,set<uint32_t>)\n";
  tax_tree.lowestCommon(tax_ids, lowest);
  for (set<uint32_t>::const_iterator t = lowest.begin(); t != lowest.end(); t++) {
    cout << "next tax_id from LCA set: " << *t << endl;
  }

  cout << "\ncalling TaxTree::lowestCommonActual(set<uint32_t>,set<uint32_t>)\n";
  lowest.clear();
  tax_tree.lowestCommonActual(tax_ids, lowest);
  for (set<uint32_t>::const_iterator t = lowest.begin(); t != lowest.end(); t++) {
    cout << "next tax_id from LCA set: " << *t << endl;
  }

  cout << endl << "paths to root for the tax_ids:\n";
  for (set<uint32_t>::const_iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
    tax_tree.printPathToRoot(*t);
  }

  cout << "\ncalling TaxTree::lowestCommon(set<uint32_t>,int)\n";
  lowest.clear();
  n = tax_tree.lowestCommon(tax_ids, 3);
  if (n != 0) {
    TaxNode *nn = const_cast<TaxNode*>(n);
    nn->write();
  }


  tax_ids.clear();
  tax_ids.insert(157976);
  tax_ids.insert(157999);
  tax_ids.insert(157975);
  tax_ids.insert(539837);
  cout << endl << "repeating test with a subset: 157976, 157999, 157975, 539837\n";
  cout << "LCA should be 984\n";

  cout << "calling TaxTree::lowestCommon(set<uint32_t>)\n";
  n = tax_tree.lowestCommon(tax_ids);
  cout << "lowestCommon returned: " << n << endl;
  if (n != 0) {
    TaxNode *nn = const_cast<TaxNode*>(n);
    nn->write();
  }

  cout << "\ncalling TaxTree::lowestCommon(set<uint32_t>,set<uint32_t>)\n";
  lowest.clear();
  tax_tree.lowestCommon(tax_ids, lowest);
  for (set<uint32_t>::const_iterator t = lowest.begin(); t != lowest.end(); t++) {
    cout << "next tax_id from LCA set: " << *t << endl;
  }

  cout << "\ncalling TaxTree::lowestCommonActual(set<uint32_t>,set<uint32_t>)\n";
  lowest.clear();
  tax_tree.lowestCommonActual(tax_ids, lowest);
  for (set<uint32_t>::const_iterator t = lowest.begin(); t != lowest.end(); t++) {
    cout << "next tax_id from LCA set: " << *t << endl;
  }

  cout << "\ncalling TaxTree::lowestCommon(set<uint32_t>,int)\n";
  lowest.clear();
  n = tax_tree.lowestCommon(tax_ids, 3);
  if (n != 0) {
    TaxNode *nn = const_cast<TaxNode*>(n);
    nn->write();
  }
}
