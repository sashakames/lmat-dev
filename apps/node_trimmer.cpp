#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "kencode.hpp"
#include "all_headers.hpp"
#include <stack>
#include <queue>
#include "TaxTree.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::set;
using std::fstream;
using std::ofstream;
using std::ostream;

using namespace kencode_ns;
using namespace metag;

static void 
usage(const string&,const string&); 

int main(int argc, char* argv[]) 
{

#if USE_GENOME_LOCATIONS != 0
#error "Must compile with USE_GENOME_LOCATIONS=0"
#endif
   char c = '\0';
   bool prn_help = false;
   string kmer_db_fn, taxtree_fn, seqids_fn, outfile, gtax_file;
   const string opt_string="t:d:t:g:h";   
   while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
      switch(c) {
      case 'g':
         gtax_file = optarg;
         break;
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      case 't':
         taxtree_fn = optarg;
         break;
         break;
      case 'h':
         prn_help = true;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
         prn_help = true;
         break;
      }
   }
   if( prn_help ) {
      usage(argv[0],opt_string);
      exit(0);
   }

   TaxTree tax_tree(taxtree_fn.c_str());
   tax_tree.findLeaves();

   TaxNodeIDArray tax_id_data;
   GenomeIdToTaxId genome_to_taxid(gtax_file.c_str());

   FILE *fp = Utils::openReadFile(kmer_db_fn.c_str());
   KmerFileMetaData metadata;
   metadata.read(fp);
   uint32_t kmer_len = metadata.kmerLength();
   uint64_t kmer_count = metadata.size();
   cout<<"KmerLen: "<<kmer_len<<" cnt: "<<kmer_count<<endl; 

   uint64_t test, sanity = ~0;

   //mark all tax nodes that have genomes associated with them;
   KmerNode kmer_node;
   KmerNode *w = &kmer_node;
   set<uint32_t> tax_ids;
   set<uint32_t> all;
   const std::vector<uint32_t> *p;
   for (uint64_t j=0; j<kmer_count; j++) {
      w->read(fp, &genome_to_taxid, &tax_id_data);  
      w->getTaxIDs(tax_ids, &tax_id_data);
      for (set<uint32_t>::const_iterator t = tax_ids.begin(); t != tax_ids.end(); t++) {
        if (all.find(*t) == all.end()) { //all is the set of all leaf nodes that
                                         //have associated genomes
          all.insert(*t);

          //march up the path to root, marking the interior nodes
          assert(tax_tree.find(*t) != tax_tree.end()); 
          tax_tree[*t]->mark = true;
          p = tax_tree[*t]->getPathToRoot();
          for (size_t j=0; j<p->size(); j++) {
            assert(tax_tree.find((*p)[j]) != tax_tree.end()); 
            tax_tree[(*p)[j]]->mark = true;
          }
        }
      }
      if ((1+j) % 1000 == 0) {
        assert(fread(&test, sizeof(uint64_t), 1, fp) == 1);
        assert(test == sanity);
      }
   }
   fclose(fp);

   TaxNode *n;
   set<uint32_t> discard_leaves;
   int discard_interior = 0;
   int keep_interior = 0;
   int leaf_count = 0;
   bool flag;
   set<uint32_t> immediate_parents;
   int parents = 0;

   //for each tax tree node ...
   for (TaxTree::const_iterator t = tax_tree.begin(); t != tax_tree.end(); t++) {
     n = t->second;
     flag = false;

     //discover if any (descendant) leaf nodes have kmers (genomes) associated with them
     for (set<uint32_t>::const_iterator t2 = n->m_leaves.begin(); t2 != n->m_leaves.end(); t2++) {
       assert(tax_tree.find(*t2) != tax_tree.end());
       if (tax_tree[*t2]->mark) {
         flag = true;
         break;
       }
     }

     //count total number of leaf nodes
     if (t->second->isLeaf()) {
       ++leaf_count;
     } 

     //deal with interior nodes
     else {
       //we can discard this node
       if (!flag) {
         ++discard_interior;
         //we can also discard all descendent leaves
         for (set<uint32_t>::const_iterator t2 = n->m_leaves.begin(); t2 != n->m_leaves.end(); t2++) {
           discard_leaves.insert(*t2);
         }

       }
     
       //must keep this node, and all descendent leaves
       else {
         ++keep_interior;
         const std::set<uint32_t>& c = t->second->getChildren();
         bool is_immediate_parent = false;
         for (set<uint32_t>::const_iterator t2 = n->m_leaves.begin(); t2 != n->m_leaves.end(); t2++) {

           //is one of the interior node's descendant leaves an immediate child?
           if (c.find(*t2) != c.end()) {
             is_immediate_parent = true;
           }
       }
       if (is_immediate_parent) {
         ++parents;
       }
     }
   }
   }

   //report results
   cout << endl;
   cout << "interior nodes: " << keep_interior+discard_interior << endl;
   cout << "interior nodes, keep: " << keep_interior << endl;
   cout << "interior nodes, discard: " << discard_interior << endl;

   cout << "leaf count: " << leaf_count << endl;
   cout << "leaf count, keep: " << (leaf_count - discard_leaves.size())<< endl;
   cout << "leaf count, discard: " << discard_leaves.size()<< endl;

   cout << "number of leaves with associated genomes: " << all.size() << endl;

   cout << "interior nodes that have children that are leaf nodes,";
   cout << "where at least one leaf node is associated with a genome: " << parents << endl;


   return 0; 
}



void
usage(const string& prog_name, const string& optstr) {
   cout<<"Usage: "<<prog_name<<" "<<optstr<<endl;
}
