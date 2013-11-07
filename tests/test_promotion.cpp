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
  // load lookup tables
  char buf[1024];
  sprintf(buf, "%s/gid2fasta_header.dat", base);
  GidToName gid_to_name(buf);
  sprintf(buf, "%s/viral_gid2tid.dat", base);
  GidToTid gi_to_tid(buf);

  // construct and populate kmer hash map
  // usr/gapps/dst/metag/data/viral_kmers_16.bin is a subset
  // contains the subset of all_virus that I was able to match
  // with names in the NCBI taxonomy (which is most of them)
  //sprintf(buf, "%s/all_virus_with_gi.fa_16.int.bin.loc", "/p/lscratchc/hysom/");
  sprintf(buf, "%s/_x.fa.bin.loc", base);
  StopWatch clock;
  clock.start();
  KmerDB kmer_hash(buf, 16);
  int K = kmer_hash.getKmerLength();
  cout << "the KmerDB contains " << kmer_hash.size()
       << " kmers of length " << K << endl;
  cout << "time to create KmerDB: " << clock.stop() << endl;

  // construct taxonomy objects
  clock.reset();
  clock.start();
  sprintf(buf, "%s/tax_nodes.dat", base);
  TaxTree tax_tree(buf);
  cout << "time to create TaxTree: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.addTidToKmerNodes(kmer_hash, gi_to_tid);
  cout << "tax tree contains " << tax_tree.size() << " nodes\n";
  cout << "time to add taxonomy node IDs to KmerNodes: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.addPathsToRoot();
  cout << "time to fill in paths from tax nodes to root: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.findChildren();
  cout << "time to find children: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.findLeaves();
  cout << "time to find leaves: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.addKmersToTaxNodes(kmer_hash);
  cout << "time to add kmers to tax nodes: " << clock.stop() << endl;

  for (std::map<int, TaxNode*>::const_iterator t = tax_tree.m_tree.begin(); t != tax_tree.m_tree.end(); t++) {
    const TaxNode* tn = t->second;
    if (tn->rank() == "species") {
      cout << tn->name() << " :: children: " << tn->m_children.size() << " :: ";
      tax_tree.printPathToRoot(tn->id());
    }
#if 0
    if (! tn->isLeaf())
      cout << tn->rank() << " " << tn->m_path_to_root.size() << " " << tn->m_path_to_root.back() << endl;
     
     if (tn->m_path_to_root.size() > 3) {
      const TaxNode *first_leaf;
      int leaf_id;
      for (set<int>::const_iterator t2 =  tn->m_leaves.begin(); t2 != tn->m_leaves.end(); t2++) {
        leaf_id = *t2;
        assert(tax_tree.m_tree.find(leaf_id) != tax_tree.m_tree.end());
        first_leaf = tax_tree.m_tree[leaf_id];
        break;
      }
      cout << "next tax node: " << tn->name() << " rank: " << tn->rank() << " distance from root: " << tn->m_path_to_root.size() << " number of leaves: " << tn->m_leaves.size() <<  " number of kmers associated with first child: " << first_leaf->m_kmers.size() << endl;

      //loop over kmers associated with first leaf"
      for (set<uint64_t>::const_iterator t3 = first_leaf->m_kmers.begin(); t3 != first_leaf->m_kmers.end(); t3++) {
        uint64_t kmer = *t3;
        int count = 0;
        //discover the number of leaves that have this kmer
        for (set<int>::const_iterator t =  tn->m_leaves.begin(); t != tn->m_leaves.end(); t++) {
          leaf_id = *t;
          assert(tax_tree.m_tree.find(leaf_id) != tax_tree.m_tree.end());
          const TaxNode* tn = tax_tree.m_tree[leaf_id];
          if (tn->m_kmers.find(kmer) != tn->m_kmers.end()) {
            ++count;
          }
        }
        cout << "  next kmer " << kmer << " is contained in " << count << " of " << tn->m_leaves.size() << " as a percentage: " << 100.0*(double)count/(double)tn->m_leaves.size() << endl;
      }
    }
#endif
    }
}

