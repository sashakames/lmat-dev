#include "../include/all_headers.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cassert>

#define MMAP_FILE "/tmp/test.mmap"
#define MMAP_SIZE ((size_t)1 << 30)

using namespace std;
using namespace metag;


char *base = "/home/hysom/metag_repo/dev/tests/";

ofstream saveme;

// returns, in kmers, the set of kmers of length K that occur
// in "read."  The kmers are canonicalized, e.g, both fwd and
// reverse complement directions are considered, and whichever
// is smaller is used
void getReadKmers(set<uint64_t> &kmers, string &read, int K);

// for each kmer, print:
//   for each genome in which it appears:
//     print the path from leaf to root in the taxonomy forest
const TaxNode * printKmerInfo(uint64_t kmers,
                   KmerDB &mer_hash,
                   GidToName &gid_to_name,
                   GidToTid &gid_to_tid,
                   TaxTree & tax_tree,
                   set<int> &was_seen);


int main(int argc, char *argv[]) {
  cerr << "starting\n";
  bool db_exists = (argc > 1);
  char *mode = (db_exists) ? "r+" : "w+";
  char buf[1024];
  KmerDB *kmer_hash = new KmerDB;

#if WITH_PJMALLOC == 1
  cout << db_exists ? "DB exists" : "DB create" << endl;
  kmer_hash->open(MMAP_FILE, mode, MMAP_SIZE);
#endif

  if (!db_exists) {
    sprintf(buf, "%s/all_virus_with_gi.fa_16.int.bin.loc", base);
    // create data base
  }

#if WITH_PJMALLOC == 1
  kmer_hash->close();
#endif

  delete kmer_hash;
  exit(0);

  // load lookup tables
/*
  sprintf(buf, "%s/gid2fasta_header.dat", base);
  GidToName gid_to_name(buf);
  sprintf(buf, "%s/viral_gid2tid.dat", base);
  GidToTid gi_to_tid(buf);
*/

  // construct and populate kmer hash map
  // usr/gapps/dst/metag/data/viral_kmers_16.bin is a subset
  // contains the subset of all_virus that I was able to match
  // with names in the NCBI taxonomy (which is most of them)
  StopWatch clock;
  clock.start();
  int K = kmer_hash->getKmerLength();
  cout << "the KmerDB contains " << kmer_hash->size()
       << " kmers of length " << K << endl;
  cout << "time to create KmerDB: " << clock.stop() << endl;

#if 0
  // construct taxonomy objects
  clock.reset();
  clock.start();
  sprintf(buf, "%s/tax_nodes.dat", base);
  TaxTree tax_tree(buf);
  cout << "time to create TaxTree: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.addTidToKmerNodes(*kmer_hash, gi_to_tid);
  cout << "tax tree contains " << tax_tree.size() << " nodes\n";
  cout << "time to add taxonomy node IDs to KmerNodes: " << clock.stop() << endl;

  clock.reset();
  clock.start();
  tax_tree.addPathsToRoot();
  cout << "time to fill in paths from tax nodes to root: " << clock.stop() << endl;

/*
  clock.reset();
  clock.start();
  tax_tree.findChildren();
  cout << "time to find children: " << clock.stop() << endl;
*/

  clock.reset();
  clock.start();
  tax_tree.findLeaves();
  cout << "time to find leaves: " << clock.stop() << endl;

  // open read file; the "reads" are subsequences of the genomes
  // used to construct the kmer DB
  //sprintf(buf, "%s/v_80_fake_reads.fa", base);
  sprintf(buf, "fake_reads_80.fa");
  ifstream in(buf);
  assert(in.is_open());

  //main loop
  string header, read;
  set<uint64_t> kmers;
  set<int> was_seen;

//saveme.open("fake_reads_80.fa");
//assert(saveme.is_open());


  while (true) {
    // get the next read; the "truth" (the genome from which the
    // read was taken) is in the header.  Of course, the read
    // might appear in other other genomes; we'll check for that later;
    // the reads do not have degenerate bases.
    header = "";
    getline(in, header);
    if (header.size() == 0) {
      break;
    }
    getline(in, read);

    cout << endl << "========================================================\n";
    cout << "header from next read: " << header << endl;
    cout << "the read: " << read << endl;

    // get the kmers in the read, then print out information about each
    getReadKmers(kmers, read, K);
    map<const TaxNode*, int> m;
    was_seen.clear();
    cout << "kmer count: " << kmers.size() << endl;
int c = 0;
    for (set<uint64_t>::const_iterator t = kmers.begin(); t != kmers.end(); t++) {
    ++c;

cout << "kmer " << c << " of " << kmers.size() << endl;
      //reads are subsequences of genomes in the viral DB,
      //so they'd better be in the table!
      if (kmer_hash->find(*t) != 0) {
      cout << "FOUND KMER\n";
        const TaxNode* tn = printKmerInfo(*t, *kmer_hash, gid_to_name, gi_to_tid, tax_tree, was_seen); 
        if (m.find(tn) == m.end()) {
          m[tn] = 0;
        }
        m[tn] += 1;
        /*
        if (printKmerInfo(*t, kmer_hash, gid_to_name, gi_to_tid, tax_tree, was_seen)) {
//        saveme << header << endl;
 //       saveme << read << endl;
        } else {
      cout << "NOT FOUND KMER\n";
      }
      */
      }
    }

    cout << endl << "all lowest common anscestors:\n";
    for (map<const TaxNode*, int>::const_iterator t = m.begin(); t != m.end(); t++) {
        cout << "# times found: " << t->second << " " << t->first->rank() << " " << t->first->name() << endl;
    }
  }
//saveme.close();
#endif
}

void getReadKmers(set<uint64_t> &kmers, string &read, int K) {
  kmers.clear();
  char work[1024];
  uint64_t kmer, rc;
  for (size_t j=0; j<read.size()-K+1; j++) {
    sprintf(work, "%s", read.substr(j, K).c_str());
    kmer = parse_dna::mer_string_to_binary(work, K);
    rc = parse_dna::reverse_complement(kmer, 16);
    if (kmer < rc) {
      kmers.insert(kmer);
    } else {
      kmers.insert(rc);
    }
  }
}

const TaxNode * printKmerInfo(uint64_t kmer,
                   KmerDB &kmer_hash,
                   GidToName &gid_to_name,
                   GidToTid &gi_to_tid,
                   TaxTree & tax_tree,
                   set<int> &was_seen) {
  KmerNode * kmer_node = kmer_hash.find(kmer);
  assert(kmer_node != 0);

  bool retval = false;

  static vector<int> nodes;


  const std::set<int> & tax_nodes = kmer_node->getTaxonomyNodes();


//  cout << "kmer is associated with " << tax_nodes.size() << " taxonomy nodes\n";
  for (set<int>::const_iterator t = tax_nodes.begin(); t != tax_nodes.end(); t++) {
    if (was_seen.find(*t) == was_seen.end()) {
      was_seen.insert(*t);
      cout << "\npath to parent from tax node: " << *t << endl;
      const TaxNode *node = tax_tree.find(*t);
      if (node == 0) {
        cout << "failed to find TaxNode for tid: " << *t << endl;
        return false;
      }

      tax_tree.getPathToRoot(*t, nodes);
      for (size_t j=0; j<nodes.size(); j++) {
        const TaxNode *nd = tax_tree.find(nodes[j]);
        if (nd == 0) {
          cout << "failed to find TaxNode for tid: " << *t << endl;
          return false;
        }
        string rank = nd->rank();
        printf("%9d %15s %30s\n", nd->id(), rank.c_str(), nd->name().c_str());
      }
    }
  }
  const set<int> *t1 = &tax_nodes;
  if (tax_nodes.size() > 1) {
    retval = true;
  }
  set<int> *tmp = const_cast<set<int>*>(t1);
  const TaxNode *lc = tax_tree.lowestCommon(*tmp);

if (lc == 0) {
  cout << "FAILED to find LCA for " << tax_nodes.size() << " TaxNodes\n";
} else {
  cout << "FOUND LCA for " << tax_nodes.size() << " taxonomy nodes" << endl;
  cout << "tid: " << lc->id() << "  rank: " << lc->rank() << "  name: " << lc->name() << endl;
  cout << "percent of LCA's leaves in the tax node set: " << 100.0*(double)tax_nodes.size()/lc->m_leaf_count << " (leaf count: " << lc->m_leaf_count << endl;
}

return lc;
//return retval;
}
