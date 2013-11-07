#ifndef __KmerNodeSorter_no_locations__
#define __KmerNodeSorter_no_locations__


#include <iostream>
#include <ostream>
#include <set>
#include <vector>
#include <map>


using namespace std;

namespace metag {

/**
 specialized KmerNodeSorter_no_locations class for sorting kmer data files
 */


class KmerNodeSorter_no_locations {
public:

  struct NodeCmp {
    bool operator() (KmerNodeSorter_no_locations* a, KmerNodeSorter_no_locations* b) { 
      return a->m_kmer < b->m_kmer;
    }  
  };

  //! write a binary representation of the kmer to file.
  void writeBinary(FILE *out) {
    assert(fwrite(&m_kmer, 8, 1, out) == 1);
    //write the number of genomes
    int c = m_genomes.size();
    assert(fwrite(&c, 4, 1, out) == 1);

    //loop over genomes 
    for (size_t j=0; j<m_genomes.size(); j++) {
      c = m_genomes[j];
      assert(fwrite(&c, 4, 1, out) == 1);
    }
  }

void  read(const char *b) {
  int g_ct, g_id; 
  uint64_t kmer;

    //read the kmer
    uint64_t *t = (uint64_t*)b;
    kmer = t[0];
    m_kmer = kmer;

    b += 8;
    int *bb = (int*)b;

    //get the number of genomes in which the kmer appears
    g_ct = bb[0];
    ++bb;

    //for each genome ...
    for (int h=0; h<g_ct; h++) {
      //get the genome ID
      g_id = bb[0];
      ++bb;
      m_genomes.push_back(g_id);
    }
}

private:

  uint64_t m_kmer;

  std::vector<int> m_genomes;

};

}

#endif
