#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <omp.h>
#include "kencode.hpp"
#include "all_headers.hpp"


#define MMAP_SIZE 0

extern int stop_early;

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;

static bool verbose=false;

size_t perm_bytes_allocd;

typedef pair<uint32_t,float> ufpair_t;
typedef pair<uint32_t,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<uint16_t,tax_data_t> label_info_t;
typedef map<uint32_t,uint32_t> hmap_t;

static char getComp(char base) {
   char c='X';
    switch (tolower(base)) {
      case 'a' : c = 't'; break;
      case 'c' : c = 'g'; break;
      case 't' : c = 'a'; break;
      case 'g' : c = 'c'; break;
      default  : c = 'X';
    }
  return c;
}

static bool isBase(char base) {
   bool res=false;
    switch (tolower(base)) {
      case 'a' : 
      case 'c' : 
      case 't' : 
      case 'g' : res=true; break;
      default  : res=false; 
    }
  return res;
}


typedef std::pair<TaxNode*,uint16_t> tncpair_t;

void updateSaveLst(const tncpair_t& chk, list<tncpair_t>& lst) {
   bool exclude=false;
   list<tncpair_t>::iterator it = lst.begin();
   const list<tncpair_t>::iterator is = lst.end();
   for(; it != is; ++it) {
      tncpair_t& curr=*it;
      if(curr.first->isAncestor(chk.first->id())) {
         //cout<<"erase "<<(*it).first->id()<<endl;
         lst.erase(it); 
      } else if( chk.first->isAncestor(curr.first->id()) ) {
         exclude=true;
      }
   }
   if(!exclude) {
      //cout<<"valid cand: "<<chk.first->id()<<endl;
      lst.push_back(chk);
   }
}

static bool isAncestor(const TaxTree& tax_tree, uint32_t prev_taxid /*ancestor */, uint32_t curr_taxid /* descendant */ ) {
   bool isOk=false;
   vector<uint32_t> path;
   TaxTree& tax_tree_tmp = const_cast<TaxTree&>(tax_tree);
   tax_tree_tmp.getPathToRoot(curr_taxid,path);
   for(unsigned p = 0; p < path.size(); ++p) {
      if(path[p] == prev_taxid) {
         isOk=true;
         break;
      }
   }
   return isOk;
}


static void updateBest(list<ufpair_t>& clst, ufpair_t conflict, const TaxTree& tax_tree, float diff_thresh) {
   while( !clst.empty() ) {
      const ufpair_t cbest = clst.front();
      if( isAncestor(tax_tree,cbest.first,conflict.first ) ) {
         //resolved conflict, we can stop
         break;
      } else if( cbest.second - conflict.second < diff_thresh ) {
        // to close a score to make a call we must remove this one
        clst.pop_front();
      } else {
        // still conflict, but too big a scoring differential to care
        break;
      }  
   }
}

static int getBest(list<ufpair_t>& cand,const vector<ufpair_t>& rank_label, const TaxTree& tax_tree, float diff_thresh ) {
   int bad_idx_start = -1;
   const int lidx = rank_label.size()-1;
   cand.push_front(rank_label[lidx]);
   for(signed i = rank_label.size()-1; i >= 0; --i) {
      if( i < lidx ) {
         const uint32_t prev_taxid = rank_label[i+1].first;
         const uint32_t curr_taxid = rank_label[i].first;
         const float  curr_taxid_score = rank_label[i].second;
         const float  prev_taxid_score = rank_label[i+1].second;
         if( !isAncestor(tax_tree,prev_taxid,curr_taxid) && !isAncestor(tax_tree,curr_taxid,prev_taxid) /* should happen but there may be a bug*/ ) {
            if(verbose) cout<<"failed here: "<<curr_taxid<<" prev "<<prev_taxid<<" "<<i<<endl;
            if( (prev_taxid_score-curr_taxid_score) < diff_thresh ) {
               bad_idx_start = i+1;
               cand.pop_front();
            } else {
               bad_idx_start = i;
            }
            break;
         } else {
            if(verbose) cout<<"Add to mix: "<<i<<" "<<rank_label[i].first<<" "<<rank_label[i].second<<endl;
            cand.push_front(rank_label[i]);
         }
      }
   }
   return bad_idx_start;
}


void
fill_in_labels(const TaxTree& taxtree, vector<uint32_t>& row, const tax_data_t& taxids, const hmap_t& idx2taxid) {
cout << "starting fill_in_labels\n";
  //   bool verbose=true;
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const uint32_t tax_val = (*mtch).second;
cout << "tax_val: " << tax_val << " count: " <<  row[tax_idx] <<endl;
      const TaxNodeHash::const_iterator row_tax_val_it = taxtree.find(tax_val); 
      if(row_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<tax_val<<endl;
         continue;
      }
      TaxNode* row_tax_val_node = (*row_tax_val_it).second;
      //cout<<"debug: "<<tax_val<<endl;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !row_tax_val_node->isLeaf() ) {
cout << "  isLeaf\n";
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         list<tncpair_t> saveLst;
         for(; it != is; ++it) {   
            const uint32_t tax_id = (*it).first;
            const uint16_t has_cnt = (*it).second;
            const TaxNodeHash::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<tax_val<<endl;
               continue;
            }
            TaxNode* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            //const TaxNode* mfill_tch = taxtree.lowestCommon(tax_id,tax_val);
            // save max value
            //cout<<"debug against: "<<tax_id<<endl;
            if( tax_val != tax_id && tax_val_node->isLeaf() && tax_val_node->isAncestor(tax_val)) {
               // want to capture the descendant tax nodes closest to row_tax_val_node 
               updateSaveLst(make_pair(tax_val_node,has_cnt),saveLst);
            }
         }
         unsigned has_cnt_sum = 0;
         list<tncpair_t>::const_iterator sit = saveLst.begin();
         const list<tncpair_t>::const_iterator sis = saveLst.end();
         for(; sit != sis; ++sit) {
            has_cnt_sum += (*sit).second;
         }
         if( has_cnt_sum > 0 ) {
            row[tax_idx] = has_cnt_sum;
cout << "    final count: " << row[tax_idx]<<endl;
            if(verbose) cout<<"Final Tax Node Val "<<tax_val<<" final transfer cnt for taxid="<<row[tax_idx]<<endl;
         }
      }
   }
   if(verbose) cout<<"fill_row: "; 
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const uint32_t tax_val = (*mtch).second;
      if(verbose) cout<<"["<<tax_val<<","<<row[tax_idx]<<"]";
   }
   if(verbose) cout<<endl;
}

struct TCmp { 
   TCmp(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<uint32_t,float>& a, const pair<uint32_t,float>& b) const {
	if( fabs(a.second- b.second) < 0.001 ) {
		const int adepth = (*_imap.find(a.first)).second;
		const int bdepth = (*_imap.find(b.first)).second;
		return adepth > bdepth;
        } 
	return a.second < b.second; }
   const hmap_t& _imap;
};

struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(0), _min_match(0) {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh;
   float _min_match;
   
};

void 
construct_labels(const TaxTree& tax_tree, const vector<label_info_t>& label_vec, const set<uint32_t>& taxid_lst, const hmap_t& tax2idx, const hmap_t& idx2taxid, ofstream& ofs, size_t mer_len, const ScoreOptions& sopt) {
cout << ">>>>>>>>>>> construct_labels\n";
   //   bool verbose = true;
   vector<vector<uint32_t> > label_matrix(label_vec.size());
   const unsigned num_tax_ids = taxid_lst.size();
   vector<bool> any_kmer_match(label_vec.size(),false) ;
   unsigned cnt_fnd_kmers=0, coverage = 0;
   int last_hit_pos = -1;
   for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
      label_matrix[pos].resize(num_tax_ids,0);
      tax_data_t::const_iterator it = label_vec[pos].second.begin();
      const tax_data_t::const_iterator is = label_vec[pos].second.end();
      unsigned any=0;
      for(; it != is; ++it) {
         const uint32_t tax_id = (*it).first;
         const uint16_t has_cnt = (*it).second;
         const unsigned idx = (*(tax2idx.find(tax_id))).second;
         label_matrix[pos][idx] = has_cnt;
         ++any;
         if(verbose) cout<<"check: "<<pos<<" "<<idx<<" "<<tax_id<<" "<<has_cnt<<endl;
      }   
      if( any > 0 ) {
         any_kmer_match[pos] = true;
cout << ">>>>>>>>>>> calling fill_in for  pos: " << pos << endl;
         fill_in_labels(tax_tree,label_matrix[pos],label_vec[pos].second,idx2taxid);
         ++cnt_fnd_kmers;

         if(last_hit_pos == -1 ) {
            coverage += mer_len;
         } else if( pos < (last_hit_pos + mer_len) ) {
            const int overlap = (last_hit_pos + mer_len) - pos;
            const int add_coverage = mer_len - overlap;
            coverage += add_coverage;
         } else { // pos >= (last_hit_pos+mer_len 
            coverage += mer_len;
         }
         last_hit_pos=pos;
      }
   }
   vector<ufpair_t> rank_label(num_tax_ids,make_pair(0,0)); 
   for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
      float tot_genome_cnt = 0, found_genome_cnt = 0;
      const uint32_t taxid = (*(idx2taxid.find(tax_idx))).second;
      if(verbose) cout<<"col: "<<tax_idx<<" "<<taxid;
      bool noMatch=false;
      for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
         if( any_kmer_match[pos] ) {
            if( sopt._strict_kmer_match && label_matrix[pos][tax_idx] == 0 ) {
               noMatch=true;
               break;
            }
            if( sopt._equal_kmer_vote ) { 
               tot_genome_cnt += 1; // equal weight;
               //const float pcnt_has = (float)label_matrix[pos][tax_idx] / (float)label_vec[pos].first;
               if( label_matrix[pos][tax_idx] > 0 ) {
               	   found_genome_cnt += 1;
               }
            } else {
               found_genome_cnt += label_matrix[pos][tax_idx];
               tot_genome_cnt += label_vec[pos].first;
            }
            if(verbose) cout<<" pos="<<pos<<" "<<label_vec[pos].first<<" "<<label_matrix[pos][tax_idx]<<" taxid="<<taxid<<" "<<found_genome_cnt<<" "<<tot_genome_cnt<<endl;
         }
      }
      float label_prob = 0.0;
      if( !noMatch ) {
         label_prob = (float)found_genome_cnt / (float)tot_genome_cnt;
      }
      if(verbose) cout<<" "<<label_prob<<endl;
      if( label_prob > sopt._min_match ) {
         rank_label[tax_idx] = make_pair(taxid,label_prob);
      } else {
         rank_label[tax_idx] = make_pair(taxid,-1);
      }
   }
   TCmp tcmp(sopt._imap);
   sort(rank_label.begin(),rank_label.end(),tcmp);
   const float pcnt_cov = (float)coverage /(float)(label_matrix.size()+mer_len);
   ofs<<"read_label "<< cnt_fnd_kmers<<" "<<pcnt_cov;
   assert(rank_label.size() > 0);  
   list<ufpair_t> valid_cand;

// get best guess and then prune back according to conflicting matches
   const int sidx = getBest(valid_cand,rank_label,tax_tree, sopt._diff_thresh);
   if(verbose) cout<<"debug: "<<sidx<<" "<<valid_cand.front().first<<" "<<valid_cand.front().second<<endl;
   for(signed i = sidx; i >= 0; --i) {
      const ufpair_t best_guess=valid_cand.front();
      const float  curr_taxid_score = rank_label[i].second;
      if( curr_taxid_score == -1 ) {
         break;
      }
      const float diff = best_guess.second - curr_taxid_score;
      if( diff <= sopt._diff_thresh ) {
         updateBest(valid_cand,rank_label[i],tax_tree, sopt._diff_thresh);
      } else {
         break;
      }
   } 
   if(sopt._prn_all) {
      for(signed i = rank_label.size()-1; i >= 0; --i) {
         ofs<<" "<<rank_label[i].first<<" "<<rank_label[i].second;
      }
   }
   if( !valid_cand.empty() ) {
      const ufpair_t best_guess = valid_cand.front();
   	ofs<<" "<<best_guess.first<<" "<<best_guess.second;
   } else {
   	ofs<<" "<<-1<<" "<<-1;
   }

   ofs<<endl;
   
}
/*
** step through each k-mer in the read in a single direction and update the "presence/absence" vector (pa_vec) when a k-mer is found
**
** h => k-mer database
** read => query read
** mer_len => k-mer length
** pa_vec => a 'bit' vector, with one element for each position in the read. Records whether the kmer for that position  is found in
**           the database.
**
** Note that the read is search in each direction separately, however, since the k-mers are stored in canonical order
** strand differentiation is not made, thus presence absence from either strand is stored in the single direction
** 
*/

void retrieve_kmer_labels(TaxTable* table, const string& read, const size_t mer_len, vector<label_info_t>& label_vec, set<uint32_t>& taxid_lst, hmap_t& tax2idx, hmap_t& idx2tax, bool revStrnd, float threshold)
{

  //  bool verbose=true;
  kencode_c ken(mer_len);
  size_t read_len = read.length();
  uint64_t kmer_id = 0x0;
  bool bad_last_kmer=true;
  unsigned next_valid_pos = 0;
  
   
  for(size_t j = 0; j < read.size(); ++j) {
    if(verbose) cout<<"retrieve kmer at posit:"<<j<<endl;
    if(!isBase(read[j])) {
        bad_last_kmer=true;
        next_valid_pos = (j + mer_len);
        continue;
    }
    if( bad_last_kmer && (j >= mer_len - 1) && j >= next_valid_pos ) {
        const string nkmer=read.substr((j-mer_len)+1,mer_len);
        kmer_id = ken.kencode(nkmer);
        if(verbose) cout<<"create new kmer at posit:"<<j<<endl;
        bad_last_kmer=false;
    } else if( !bad_last_kmer) {
        kmer_id = ken.kencode(read[j]);
        bad_last_kmer=false;
        if(verbose) cout<<"update kmer at posit:"<<j<<endl;
    }
    size_t pos = j - mer_len + 1;
    if (revStrnd) {
       pos = read_len - j - 1;
    }
    if(!bad_last_kmer && label_vec[pos].first == 0) {
       TaxNodeStat *h = new TaxNodeStat(*table);
       if(verbose) cout<<"lookup kmer at posit:"<<j<<endl;
       h->begin(kmer_id);
       unsigned dcnt = 0, mtch = 0;
       while( h->next() ) {  
         const uint32_t tid = h->taxid();
         const uint16_t ng = h->taxidCount();

         if(dcnt==0) {
            assert( ng > 0 );
            label_vec[pos].first = ng;
         }
         //collect the unique set of tax ids
         const uint16_t pr_cnt = h->present();
         const float pcnt = (float)pr_cnt / (float) label_vec[pos].first;
         if( pcnt >= threshold ) {
            const uint16_t pr_cnt = h->present();
            //const float pcnt = (float)pr_cnt / (float) label_vec[pos].first;
            label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
            if( tax2idx.find(tid) == tax2idx.end() ) {
               const unsigned idx = taxid_lst.size();
               tax2idx[tid] = idx;
               idx2tax[idx] = tid;
               taxid_lst.insert(tid);
            }
         // Note: ng should not change with multiple calls to next here
            if(verbose) {
               if( dcnt == 0 ) cout<<"debug kmer "<<kmer_id<<" gc="<<ng;
               cout<<" ["<<pos<<" "<<tid<<" "<<pr_cnt<<"]" ;
            }
            dcnt++;
         } else {
            //if(verbose && pcnt >= 0.05) cout<<" [failed "<<pos<<" "<<tid<<" "<<pr_cnt<<"]" ;
         }
         ++mtch; 
       }
       //}
       if( verbose ) cout<<"num taxids: "<<taxid_lst.size();
       if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
       if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
       delete h;
    }
  }
}

void proc_line(const TaxTree& tax_tree, int ri_len, string &line, int k_size, TaxTable *table, ofstream &ofs, float threshold, const ScoreOptions& sopt) {

  int thread_num = omp_get_thread_num();
  assert(ri_len >= 0 && ri_len <= line.length());

    for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
       string rval;
       if( ri_len < line.length() ) {
         rval = line.substr(chk,ri_len);
       } else {
         rval = line;
       }
       const unsigned read_len = rval.length();
       if( static_cast<signed>(read_len) < ri_len ) {
         break;
       }  
       if( read_len < k_size ) {
         if( chk == 0 ) /// read too short
            ofs<<"read_label 0 -1"<<endl; 
         break;
       } 
       vector<label_info_t> label_vec(read_len-k_size+1);
       set<uint32_t> taxid_lst; 
       hmap_t tax2idx, idx2tax;

cout << "calling retrieve_kmer_labels: "<< rval << endl;
       retrieve_kmer_labels(table, rval, k_size, label_vec, taxid_lst, tax2idx, idx2tax, false, threshold);

cout << "taxid_lst size: " << taxid_lst.size() << endl;
cout << "idx2tax size: " << idx2tax.size() << endl;

       string rev_cmp(rval.size(),'\0');
       assert(rval.size()>0);
       for (int j=(rval.size()-1); j >= 0; --j) {
         const int bidx = (rval.size()- 1) - j;
         rev_cmp[j] = getComp(rval[bidx]);
       }
       retrieve_kmer_labels(table, rev_cmp, k_size, label_vec, taxid_lst, tax2idx, idx2tax, true, threshold);

cout << "taxid_lst size: " << taxid_lst.size() << endl;
cout << "idx2tax size: " << idx2tax.size() << endl;

       //horribly slow probably
       //find LCA for retrieved nodes to make sure it's part of the set
       if( !taxid_lst.empty() ) {

         if( taxid_lst.size() > 1 ) {
            const TaxNode* lca = tax_tree.lowestCommon(taxid_lst, thread_num);
            if(!lca) {
               if(verbose) {
                  cerr<<"Error did not get lca for: "<<line<<endl;
                  /*
                  it = taxid_lst.begin();
                  for(; it != is; ++it) {
                     cerr<<"taxid="<<*it<<endl;
                  }
                  */
               }
            } else {
               const int lca_taxid = lca->id();
cout << "adding lca: " << lca->id() << endl;
               if( taxid_lst.find(lca_taxid) == taxid_lst.end() ) {
                  assert( tax2idx.find(lca_taxid) == tax2idx.end() );
                  const unsigned idx = taxid_lst.size();
                  tax2idx[lca_taxid] = idx;
                  idx2tax[idx] = lca_taxid;
                  taxid_lst.insert(lca_taxid);
               }
            }
         }
         construct_labels(tax_tree,label_vec,taxid_lst,tax2idx,idx2tax,ofs,k_size,sopt); 
      } else {
         ofs<<"read_label 0 -1"<<endl; 
      }
    }
}


size_t *split_file(int n_threads, ifstream &file)
{
  size_t * arr = new size_t[n_threads+1];

  arr[0] = 0;

  file.seekg(0, ios::end);
  size_t end = file.tellg();

  for (size_t i = 1; i<n_threads; i++)  {

    file.seekg( i * (end / n_threads ));
    
    string junk;
    getline(file, junk);
    
    arr[i] = file.tellg();
  }

  arr[n_threads] = end;

  return arr;

}

void usage(char *execname)
{
  cout << "Usage:\n" ;
  cout << execname << " -d <input db file (list)> -i <query fasta file> -t <number of threads> -o <output path> [-l read length]\n";
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-a:ascii input format]\n";
  cout << "note: -r makes -k unnecessary\n"; 
}


int main(int argc, char* argv[]) 
{
   char c = '\0';
   int k_size=-1, ri_len = -1;
   int n_threads = 0;

   bool restore = false;


   bool ascii = false;

   float threshold = 0.0;


   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, gid_to_tid_fn, depth_file;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   char * dumpfile = NULL;


   size_t mmap_size = 0;
   
   int max_reads =0;

   while ((c = getopt(argc, argv, "h:j:b:ye:wmpk:c:v:k:i:d:l:t:s:r o:x:f:g:z:a ")) != -1) {
      switch(c) {

      case 'h':
         stop_early = atoi(optarg);
         break;
      case 's':
         mmap_size = atoi(optarg);
         mmap_size = mmap_size * (1<<30);
         cout << "Input heap size: " << mmap_size << endl;
         break;
      case 'j':
         sopt._min_match = atof(optarg);
         break;
      case 'b':
         sopt._diff_thresh = atof(optarg);
         break;
      case 'y':
         verbose = true;
         break;
      case 'e':
         depth_file = optarg;
         break;
      case 'q':
	max_reads = atoi(optarg);
	break;
      case 'w':
         sopt._equal_kmer_vote = true;
         break;
      case 'm':
         sopt._strict_kmer_match = true;
         break;
      case 'p':
         sopt._prn_all = true;
         break;   
      case 'z':
	verbose = true;
	break;
      case 'r':
        restore = true;
        break;
      case 'a':
	ascii = true;
	break;
      case 't':
        n_threads = atoi(optarg);
        omp_set_num_threads(n_threads);
        break;
      case 'v':
         threshold = atof(optarg);
         break;
      case 'l':
         ri_len = atoi(optarg);
         break;
      case 'c':
         tax_tree_fn = optarg;
	 break;
      case 'k':
         k_size = atoi(optarg);
         break;
      case 'g':
         gid_to_tid_fn = optarg;
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      case 'o':
         ofbase = optarg;
         break;
      case 'x':
	dumpfile = optarg;
	break;
	//      case 'f':
	//	tax_tree_options = optarg;
	//	break;

      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

   if (depth_file == "" ||  ofbase == "" || n_threads == 0 || kmer_db_fn == "" || query_fn == "")  {
     cout<<ofbase<<" "<<n_threads<<" "<<kmer_db_fn<<" "<<query_fn<<" "<<depth_file<<endl; 
     usage(argv[0]);
     return -1;

   }
   if (!restore && k_size == -1)  {
     cout << "missing kmer size!\n";
     usage(argv[0]);
     return -1;

   }
   cout << "Start kmer DB load\n";
   TaxTable *taxtable;


#if USE_BOOST == 1
     bip::managed_mapped_file mfile(bip::open_only, kmer_db_fn.c_str());
     std::pair<TaxTable*, std::size_t> ret = mfile.find<TaxTable>("KmerDB");

     taxtable = ret.first;
     k_size = taxtable->get_kmer_length();
#else


#if WITH_PJMALLOC == 1

   if (restore) {

     perm(&taxtable, sizeof(taxtable));
     mopen(kmer_db_fn.c_str(), "r", mmap_size);

     if (k_size < 1)
       k_size = taxtable->get_kmer_length();

   } else 
    
#endif

{


     taxtable = new TaxTable;

     ifstream qifs(kmer_db_fn.c_str());
     if( !qifs ) {
       cerr<<"Unable to open: "<<kmer_db_fn<<endl;
       return -1;


     }



     string fname;
     
     while(qifs>>fname) {
       cout<<"register file: "<<fname<<endl;
       taxtable->registerFile(fname.c_str());
     }
     taxtable->ingest(ascii);


   }

#endif
   //cout << "End kmer DB load\n";
   //cout << "DB size is " << table->size() << endl;


   assert(k_size > 0 );
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */

   size_t *arr = split_file(n_threads, ifs);
   ifs.close();

   bool finished;

   size_t pos = 0 ;

   ofstream ofs;





   /*
   if (tax_tree_options.size()) {
    cout << "info: loading tax tree options\n";
    tax_tree.loadLcaOptions(tax_tree_options.c_str());
    //tax_tree.printLcaOptions();
    cout << "info: DONE loading tax tree options\n";
  }
   

   if (dumpfile)
     {
       tax_tree.dump(dumpfile);
       cout<<"Taxonomy tree dump finished." << endl ;
       exit(0);
     }
   */

   cout<<"Read taxonomy tree: "<<tax_tree_fn<<endl;
   TaxTree tax_tree(tax_tree_fn.c_str(), gid_to_tid_fn.c_str(), n_threads);
#if 0
   set<uint32_t> taxset;
   taxset.insert(974);
   taxset.insert(13910);
   taxset.insert(48497);
   taxset.insert(157975);
   taxset.insert(157976);
   taxset.insert(157999);
   taxset.insert(190173);
   taxset.insert(273213);
   taxset.insert(295518);
   taxset.insert(539837);
   taxset.insert(697273);
   int thread_num=0;
   const TaxNode* lca = tax_tree.lowestCommon(taxset, thread_num);
#endif

   cout<<"Read taxonomy depth: "<<depth_file<<endl;
   ifstream ifs1(depth_file.c_str());
   if(!ifs1) {
	cerr<<"unable to open: "<<depth_file<<endl;
	return -1;
   }
   uint32_t taxid,depth;
   while(ifs1>>taxid>>depth) {
	sopt._imap[taxid] = depth;
   }
   

   StopWatch clock;
   clock.start();
  
   size_t read_count = 0; 

#pragma omp parallel shared(arr, k_size, query_fn, ofbase, taxtable, tax_tree, max_reads, sopt, ri_len)  private(ifs, finished, pos, ofs, ofname, line, read_count)

   {

     read_count = 0;

     

     finished = false;
     int use_len = -1; 

     cout<<"Read query file: "<<query_fn<<endl;
     ifs.open(query_fn.c_str());
     if(!ifs) {
         cerr<<"did not open for reading: "<<query_fn<<endl;
         exit(-1);
     }
     
     ofname = ofbase;

     std::stringstream outs;

     outs << omp_get_thread_num();
     ofname += outs.str();

     ofname += ".out" ;

     ofs.open(ofname.c_str());

     ifs.seekg(arr[omp_get_thread_num()]);

     string read_buff;

     while (!finished)   {

       getline(ifs, line);

       pos = ifs.tellg();

       if ((pos >= arr[1+omp_get_thread_num()]) || (pos == -1)) {
	      finished = true;
       } 

       //#pragma omp critical    
       //cout << omp_get_thread_num()  << " - pos:" << ifs.tellg() << "\n"; 

       if (line[0] != '>' && line.length() > 1) {
          read_buff += line;
       } 
       if( (line[0] == '>' || finished) && read_buff.length() > 0 ) {
          if( ri_len == -1 ) {
            // if not specified use first read as the fixed length
            use_len = read_buff.length();
          } else {
            use_len = ri_len;
          } 
	       //#pragma omp critical
          ofs<<"read: "<<read_buff<< " pos: " << ifs.tellg() <<  endl;
       //split read into multiples of "ri_len" (read interval length)
	  proc_line(tax_tree, use_len, read_buff, k_size, taxtable, ofs, threshold,sopt);
          read_buff="";
	  read_count ++;
	  if (read_count == max_reads)
	   finished = true;

       }
     }
     ofs.close();
   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
