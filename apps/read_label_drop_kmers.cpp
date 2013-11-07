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

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::istringstream;

using namespace kencode_ns;
using namespace metag;

static bool verbose=false;

bool add_root_on_kmer_drop = true;

size_t perm_bytes_allocd;

typedef pair<uint32_t,float> ufpair_t;
typedef pair<uint32_t,uint16_t> tax_elem_t;
typedef set<tax_elem_t> tax_data_t;
typedef pair<uint16_t,tax_data_t> label_info_t;
typedef map<uint32_t,uint32_t> hmap_t;
typedef map<uint32_t,float> ufmap_t;
typedef map<uint16_t, ufmap_t> u_ufmap_t;

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

static int getReadLen(int rl) {
   int rval=0;
   if(rl < 80 ) {
     rval=50;
   } else if( rl < 100) {
      rval=80;
   } else {
      rval = 100;
   }
   return rval;
}

typedef std::pair<TaxNode*,uint16_t> tncpair_t;

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


struct CmpDepth {
   CmpDepth(const hmap_t& imap) : _imap(imap) {}
   bool operator()(const pair<uint32_t,float>& a, const pair<uint32_t,float>& b) const {
      const int adepth = (*_imap.find(a.first)).second;
      const int bdepth = (*_imap.find(b.first)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};

struct CmpDepth1 {
   CmpDepth1(const hmap_t& imap) : _imap(imap) {}
   bool operator()(uint32_t a, uint32_t b) const {
      const int adepth = (*_imap.find(a)).second;
      const int bdepth = (*_imap.find(b)).second;
      return adepth > bdepth;
   }
   const hmap_t& _imap;
};



enum match_t { eDirectMatch, eMultiMatch, ePartialMultiMatch, eNoMatch, eNoLCAError };
static string match2str(match_t match) {
   string str = "error";
   switch(match) {
   case eDirectMatch : 
      str="DirectMatch";
      break;
   case eMultiMatch : 
      str="MultiMatch";
      break;
   case ePartialMultiMatch : 
      str="PartialMultiMatch";
      break;
   case eNoMatch : 
      str="NoMatch";
      break;
   case eNoLCAError : 
      str="LCA_ERROR";
      break;
   //defualt:
      //cerr<<"Unexpected error: "<<match<<endl;
      //break;
   }
   return str;
}

static bool addToCandLineage(const ufpair_t cand, list<ufpair_t>& lineage, const hmap_t& dmap, const TaxTree& tax_tree) {
   bool addNode = false;
   if( lineage.empty() ) {
      addNode=true;
   } else {
      const hmap_t::const_iterator mtch0 = dmap.find(cand.first);
      assert(mtch0 != dmap.end());
      const unsigned cand_depth = (*mtch0).second;
      addNode = true;
      list<ufpair_t>::iterator it = lineage.begin();
      const list<ufpair_t>::iterator is = lineage.end();
      for(; it != is; ++it) {
         const uint32_t taxid = (*it).first;
         const hmap_t::const_iterator mtch = dmap.find(taxid);
         assert(mtch != dmap.end());
         const unsigned chk_depth = (*mtch).second;
         if( chk_depth > cand_depth && !isAncestor(tax_tree,cand.first,taxid) ) {
            addNode=false;
            break;
         } else if( chk_depth < cand_depth && !isAncestor(tax_tree,taxid,cand.first)) {
            addNode=false;
            break;
         } else if( chk_depth == cand_depth ) {
            addNode=false;
            break;
         }
      }
   }
   if( addNode ) {
      if(verbose) cout<<"Add to Lineage: "<<cand.first<<" "<<cand.second<<endl;
      lineage.push_back(cand);
   }
   return addNode; 
}

static bool cmpCompLineage(ufpair_t cand, const vector<ufpair_t>& lineage, set<uint32_t>& no_good, float diff_thresh, const TaxTree& tax_tree) {
   const float undef = -10000;
   bool keep_going=true;
   for(unsigned i = 0; i < lineage.size(); ++i) {
      if( isAncestor(tax_tree, lineage[i].first, cand.first )) {
         break;
      } 
      if( lineage[i].second != undef &&  (lineage[i].second - cand.second) >= diff_thresh  ) {
         if(verbose) cout<<"competing lineage is too far to care: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         keep_going=false;
         break;
      }
      if( (lineage[i].second - cand.second) < diff_thresh ) {
         if(verbose) cout<<"competing lineage is too close: "<<cand.first<<" "<<cand.second<<" "<<lineage[i].first<<" "<<lineage[i].second<<" "<<diff_thresh<<endl;
         no_good.insert(lineage[i].first);
      }
   }
   return keep_going;
}

static pair<ufpair_t,match_t> findReadLabelVer2(const vector<ufpair_t>& rank_label, float diff_thresh, const TaxTree& tax_tree, const hmap_t& taxid2idx, list<ufpair_t>& cand_lin, const hmap_t& dmap, const ufmap_t& all_cand_set) {
   match_t match = eNoMatch;
   unsigned lowest_depth = 0, highest_depth = 0;
   ufpair_t lowest=make_pair(0,0), highest = make_pair(0,0);
   signed lidx = -1;
   assert( cand_lin.empty() );
   for(signed i = rank_label.size()-1; i >= 0; --i) {
      if( !addToCandLineage(rank_label[i],cand_lin,dmap, tax_tree) )  {
         lidx = i;
         break;
      } else {
         const hmap_t::const_iterator mtch = dmap.find(rank_label[i].first);
         if( (*mtch).second> lowest_depth ) {
            lowest = rank_label[i];
            lowest_depth = (*mtch).second;
         }
         if( (*mtch).second < highest_depth || i == (signed)rank_label.size()-1 ) {
            highest = rank_label[i];
            highest_depth = (*mtch).second;
         }
      }
   }
   set<uint32_t> add_set;
   if( highest_depth != 0 ) {
      vector<uint32_t> path;
      TaxTree& tax_tree_tmp = const_cast<TaxTree&>(tax_tree);
      tax_tree_tmp.getPathToRoot(highest.first,path);
      for(unsigned i = 0; i < path.size(); ++i) {
         add_set.insert( path[i] );
         ufmap_t::const_iterator mtch = all_cand_set.find( path[i] );
         if( mtch != all_cand_set.end() ) {
            ufpair_t val = make_pair( path[i], (*mtch).second);
            cand_lin.push_back(val);
         } else {
            //compute log_odds
            const float undef = -10000;
            cand_lin.push_back( make_pair(path[i], undef) );
         }
      }
      // fill in lineage up to root here
   } 
   vector<ufpair_t> cand_lin_vec(cand_lin.size());
   list<ufpair_t>::const_iterator it = cand_lin.begin();
   const list<ufpair_t>::const_iterator is = cand_lin.end();
   for(unsigned i= 0; i < cand_lin_vec.size(); ++it, ++i) {
      cand_lin_vec[i] = *it;
   }
   CmpDepth cd(dmap);
   sort(cand_lin_vec.begin(),cand_lin_vec.end(),cd);

   set<uint32_t> no_good;

   for(signed i = lidx; i >= 0; --i) {
      if( add_set.find( rank_label[i].first ) == add_set.end() ) {
         if( !cmpCompLineage(rank_label[i],cand_lin_vec,no_good,diff_thresh,tax_tree) )  {
            if(verbose) cout<<" competing taxid too low scoring, quit search "<<rank_label[i].first<<" "<<rank_label[i].second<<endl;
            break;
         }
      }
   }
   ufpair_t taxid_call;
   if( cand_lin.empty() && no_good.empty() ) {
      match = eNoMatch;
   } else if( !cand_lin.empty() && no_good.empty() ) {
      taxid_call = lowest; 
      match = eDirectMatch;
   } else {
      vector<ufpair_t> cand_vec(cand_lin.size());
      list<ufpair_t>::const_iterator it = cand_lin.begin(); 
      const list<ufpair_t>::const_iterator is = cand_lin.end(); 
      unsigned cnt= 0;
      for(cnt= 0; it != is; ++it, ++cnt) {
         if(verbose) cout<<"merging cand_lst "<<cnt<<" "<<(*it).first<<" "<<(*it).second<<endl;
         cand_vec[cnt] = *it;
      }  
      CmpDepth cd(dmap);
      sort(cand_vec.begin(),cand_vec.end(),cd);
      float sval = 0, max_val = -10000;
      pair<uint32_t,bool> res = make_pair(0,false);
      int root_idx = -1;
      for(unsigned i = 0; i < cand_vec.size(); ++i) {
         const uint32_t tax_i = cand_vec[i].first;
         max_val = std::max(cand_vec[i].second,max_val);
         sval += cand_vec[i].second;
         ++cnt; 
         if( no_good.find(tax_i) == no_good.end() ) {
            res = make_pair( cand_vec[i].first, true );
            root_idx = i;
            break;
         }
      }
      if(!res.second) {
         taxid_call = make_pair(0,-1);
         match = eNoLCAError; 
      }  else {
         const uint32_t lca_tid = res.first;
         match = eMultiMatch; 
         if( all_cand_set.find(lca_tid) != all_cand_set.end() ) {
            assert(root_idx != -1);
            if( max_val < cand_vec[root_idx].second ) {
               match = ePartialMultiMatch;
               max_val = cand_vec[root_idx].second;
            }
         }
         taxid_call = make_pair(lca_tid, max_val);
      }
   }
   return make_pair(taxid_call,match);
}

void
fill_in_labels(const TaxTree& taxtree, vector<uint32_t>& row, const tax_data_t& taxids, const hmap_t& idx2taxid, const hmap_t& taxid2idx) {
   for(unsigned tax_idx = 0; tax_idx < row.size(); ++tax_idx) {
      const hmap_t::const_iterator mtch = idx2taxid.find(tax_idx);
      assert(mtch != idx2taxid.end());
      const uint32_t col_tax_val = (*mtch).second;
      const TaxNodeHash::const_iterator col_tax_val_it = taxtree.find(col_tax_val); 
      if(col_tax_val_it == taxtree.end()) {
         cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
         continue;
      }
      TaxNode* col_tax_val_node = (*col_tax_val_it).second;
      // see if this id's genome count needs to be update
      //  based on its children   - can skip if its a leaf
      if( !col_tax_val_node->isLeaf() && row[tax_idx] == 0 ) {
         tax_data_t::const_iterator it = taxids.begin();
         const tax_data_t::const_iterator is = taxids.end();
         unsigned has_cnt_sum = 0;
         for(; it != is; ++it) {   
            const uint32_t tax_id = (*it).first;
            const TaxNodeHash::const_iterator tax_val_it = taxtree.find(tax_id); 
            if(tax_val_it == taxtree.end()) {
               cout<<"we have a taxtree mapping problem: "<<col_tax_val<<endl;
               continue;
            }
            TaxNode* tax_val_node = (*tax_val_it).second;
            // test if tax_id is a descendant of tax_val
            // save max value
            const hmap_t::const_iterator idx_mtch = taxid2idx.find(tax_id);
            assert(idx_mtch != taxid2idx.end());
            const uint32_t tax_val_idx = (*idx_mtch).second;
            
            if( col_tax_val != tax_id && row[tax_val_idx] > 0 && tax_val_node->isAncestor(col_tax_val)) {
               // want to capture the descendant tax nodes closest to col_tax_val_node 
               has_cnt_sum = 1;
               break;
            }
         }
         if( has_cnt_sum > 0 ) {
            row[tax_idx] = has_cnt_sum;
            if(verbose) cout<<"Final Tax Node Val "<<col_tax_val<<" final transfer cnt for taxid="<<row[tax_idx]<<endl;
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
		return adepth < bdepth;
        } 
	return a.second < b.second; }
   const hmap_t& _imap;
};

struct ScoreOptions {
   ScoreOptions(hmap_t& imap) : _equal_kmer_vote(false), _strict_kmer_match(false), _prn_all(false), _imap(imap), _diff_thresh(0) /*, _min_match(0)*/ {}
   bool _equal_kmer_vote;
   bool _strict_kmer_match;
   bool _prn_all;
   hmap_t& _imap;
   float _diff_thresh;
   u_ufmap_t _rand_hits; 
};

static void loadRandHits(const string& file_lst, u_ufmap_t& rand_hits_all) {
   ifstream ifs_lst(file_lst.c_str());
   if( !ifs_lst) {
      cerr<<"Unexpected reading error: "<<file_lst<<endl;
      return;
   }
   int read_len;
   string file;
   while(ifs_lst>>read_len>>file) {
      cout<<"load: "<<read_len<<" "<<file<<endl;
      ifstream ifs(file.c_str());
      if( !ifs) {
         cerr<<"Unexpected reading error: "<<file<<endl;
         continue;
      }
      //while(ifs>>taxid>>cutoff>>dum1>>dum2>>dum3) {
      const int buff_size=20004;
      char buff[buff_size];
      while(ifs.getline(buff,buff_size)) {
         istringstream istrm(buff);
         int32_t taxid;
         float min_val=0,max_val=0,mean_val=0,stdev=0;
         istrm>>taxid>>min_val>>max_val>>mean_val>>stdev;
         //const float cutoff = mean_val+(3.0*stdev); // can play around with how cutoff is set
         const float cutoff = max_val;
         rand_hits_all[read_len][taxid]=cutoff;
      }
   } 

}

static float log_odds_score(float label_prob, float random_prob) {
      if( random_prob > 1 ) random_prob = 1;
      const float mod_prob = (1-label_prob) <= 0 ? 0.0001 : 1-label_prob;
      const float denom = mod_prob*random_prob; // must be > zero
      const float numer =  label_prob*(1-random_prob);
      if(verbose)  cout<<"log_odds_calc: "<<numer<<" "<<denom<<endl;
      const float log_odds = log( numer / denom );
      return log_odds;
}

void 
construct_labels(const TaxTree& tax_tree, const vector<label_info_t>& label_vec, const list<uint32_t>& taxid_lst, const hmap_t& tax2idx, const hmap_t& idx2taxid, ofstream& ofs, size_t mer_len, const ScoreOptions& sopt, const vector<bool> & is_forward) {
   const unsigned num_tax_ids = taxid_lst.size();
   vector<bool> any_kmer_match(label_vec.size(),false) ;
   unsigned cnt_fnd_kmers=0, coverage = 0;
   int last_hit_pos = -1;
   vector<vector<uint32_t> > label_matrix(label_vec.size());
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
      any_kmer_match[pos] = !label_vec[pos].second.empty();
      if( any_kmer_match[pos] ) {
            //fill_in_labels(tax_tree,label_matrix[pos],label_vec[pos].second,idx2taxid,tax2idx);
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
      const unsigned read_len = label_matrix.size()+ mer_len - 1;
      const  int read_len_match = getReadLen(read_len);
      const float pcnt_cov = (float)coverage /(float)(read_len);
      vector<ufpair_t> rank_label(num_tax_ids,make_pair(0,0)); 
      unsigned sig_hits = 0;
      float log_sum=0.0;
      ufmap_t all_cand_set;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         float tot_genome_cnt = 0, found_genome_cnt = 0;
         const uint32_t taxid = (*(idx2taxid.find(tax_idx))).second;
         if(verbose) cout<<"col: "<<tax_idx<<" "<<taxid;
         bool noMatch=false;

         unsigned add_coverage = 0;
         signed tax_last_hit_pos = -1;
         for(unsigned pos = 0; pos < label_vec.size(); ++pos) {
            if( any_kmer_match[pos] ) {
               if( sopt._strict_kmer_match && label_matrix[pos][tax_idx] == 0 ) {
                  noMatch=true;
                  break;
               }
               tot_genome_cnt += 1; // equal weight;
               if( label_matrix[pos][tax_idx] > 0 ) {
                     found_genome_cnt += 1;
               }
               int local_add=0;
               if( label_matrix[pos][tax_idx] > 0 ) {
                  int overlap = -1;
                  if(tax_last_hit_pos == -1 ) {
                     local_add = mer_len;
                  } else if( pos < (tax_last_hit_pos + mer_len) ) {
                     overlap = (tax_last_hit_pos + mer_len) - pos;
                     if((int)mer_len >= overlap) {
                        local_add = (mer_len - overlap);
                     }
                  } else { // pos >= (last_hit_pos+mer_len 
                     local_add = mer_len;
                  }
                  add_coverage += local_add;
                  //cout<<"ddebug: "<<local_add<<" "<<add_coverage<<" "<<overlap<<" "<<mer_len<<" "<<tax_last_hit_pos<<endl;
                  tax_last_hit_pos=pos;
               }
               if(verbose) cout<<" pos="<<pos<<" "<<label_vec[pos].first<<" "<<label_matrix[pos][tax_idx]<<" taxid="<<taxid<<" "<<found_genome_cnt<<" "<<tot_genome_cnt<<" "<<local_add<<" "<<add_coverage<<endl;
            }
         }
         float label_prob = 0.0;
         if( !noMatch ) {
            //label_prob = (float)found_genome_cnt / (float)tot_genome_cnt;
            if( sopt._equal_kmer_vote) {
               label_prob = (float)found_genome_cnt / (float) label_matrix.size();
            } else {
               label_prob = (float)add_coverage / (float)read_len;
            }
            //cout<<"debug: "<<add_coverage<<" "<<read_len<<" "<<found_genome_cnt<<" "<<tot_genome_cnt<<endl; 
         }
         if(verbose) cout<<" "<<label_prob<<endl;
         u_ufmap_t::const_iterator mtch = sopt._rand_hits.find(read_len_match);
         float random_prob = 0.0001;
         if( mtch != sopt._rand_hits.end() ) {
            const ufmap_t& rand_hits = (*mtch).second;
            if( rand_hits.find(taxid) != rand_hits.end()) {
               random_prob = (*(rand_hits.find(taxid))).second+0.0001;
            }
         }
         const float log_odds = log_odds_score(label_prob,random_prob);
         rank_label[tax_idx] = make_pair(taxid,log_odds);
         all_cand_set.insert( rank_label[tax_idx] );
         log_sum += log_odds;
         sig_hits++;
         if(verbose) cout<<"Sig match chec: "<<taxid<<" "<<rank_label[tax_idx].second<<" "<<log_odds<<" "<<label_prob<<" "<<random_prob<<endl;
      }
      list<ufpair_t> valid_cand;
      pair<ufpair_t,match_t> res = make_pair(make_pair(0,0),eNoMatch);   
      string matchType=match2str(res.second);
      const float log_avg = sig_hits > 0 ? log_sum / (float)sig_hits : 0;
      float log_std=0, max_score = 0;
      bool first=true;
      for(unsigned tax_idx = 0; tax_idx < num_tax_ids; ++tax_idx) {
         if( rank_label[tax_idx].second > 0 ) { 
            if( rank_label[tax_idx].second > max_score || first ) {
               max_score = rank_label[tax_idx].second;
               first=false;
            }
            const float val = log_avg - rank_label[tax_idx].second;
            log_std += (val*val);
         }
      }
      const float stdev = sig_hits > 1 ? sqrt(log_std/(sig_hits-1)) : 0;
      if( sig_hits > 0 ) {
         TCmp tcmp(sopt._imap);
         sort(rank_label.begin(),rank_label.end(),tcmp);
         ofs<<"read_label "<< cnt_fnd_kmers<<" "<<pcnt_cov<<" "<<max_score<<" "<<log_avg<<" "<<stdev<<" ";
         res = findReadLabelVer2(rank_label,stdev,tax_tree,tax2idx,valid_cand,sopt._imap,all_cand_set);    
         if(sopt._prn_all) {
            for(signed i = rank_label.size()-1; i >= 0; --i) {
               if( rank_label[i].second >= 0  || verbose ) {
                  ofs<<" "<<rank_label[i].first<<" "<<rank_label[i].second;
               }
            }
         }
         matchType=match2str(res.second);
      }
      if( res.second == eDirectMatch) {
         const ufpair_t best_guess = res.first;
         ofs<<" "<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }  else if( res.second == eMultiMatch || res.second == ePartialMultiMatch ) {
         if( !sopt._prn_all) {
            list<ufpair_t>::const_iterator it = valid_cand.begin();
            const list<ufpair_t>::const_iterator is = valid_cand.end();
            for( ; it != is; ++it) {
               ofs<<" "<<(*it).first<<" "<<(*it).second;
            }
         }
         const ufpair_t best_guess = res.first;
         ofs<<" "<<best_guess.first<<" "<<best_guess.second<<" "<<matchType;
      }
      else if( eNoMatch ) {
         ofs<<" "<<-1<<" "<<-1<<" "<<matchType;
      } else {
         cerr<<"Unexpected match type"<<endl;
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

   void retrieve_kmer_labels(TaxTable* table, const string& read, const size_t mer_len, vector<label_info_t>& label_vec, list<uint32_t>& taxid_lst, hmap_t& tax2idx, hmap_t& idx2tax, bool revStrnd, float threshold, vector<bool>& is_forward, const hmap_t& dmap, const TaxTree& tax_tree, uint16_t max_count) 
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
          list<uint32_t> obs_tids;
          while( h->next() ) {
            uint32_t tid = h->taxid();
            if( tid == 20999999) continue;
            const uint16_t ng = h->taxidCount();


            //start: hack to drop kmers with too many tids
            if (ng > max_count) {
              if (add_root_on_kmer_drop ) { //start: assign root
                tid = 1;
                label_vec[pos].first = tid; 

                const uint16_t pr_cnt = tid; //please fixme
                obs_tids.push_back(tid);
                label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
                is_forward[pos] = !revStrnd;
                if( tax2idx.find(tid) == tax2idx.end() ) {
                  const unsigned idx = taxid_lst.size();
                  tax2idx[tid] = idx;
                  idx2tax[idx] = tid;
                  taxid_lst.push_back(tid);
                }
              } //end: assign root
              break;
            }  //end: hack to drop kmers with too many tids

            if(dcnt==0) {
               assert( ng > 0 );
               label_vec[pos].first = ng;
            }
            //collect the unique set of tax ids
            const uint16_t pr_cnt = h->present();
               obs_tids.push_back(tid);
               label_vec[pos].second.insert( make_pair(tid,pr_cnt) );
               is_forward[pos] = !revStrnd;
               if( tax2idx.find(tid) == tax2idx.end() ) {
                  const unsigned idx = taxid_lst.size();
                  tax2idx[tid] = idx;
                  idx2tax[idx] = tid;
                  taxid_lst.push_back(tid);
               }
            // Note: ng should not change with multiple calls to next here
               if(verbose) {
                  if( dcnt == 0 ) cout<<"debug kmer "<<kmer_id<<" gc="<<ng;
                  cout<<" ["<<pos<<" "<<tid<<" "<<pr_cnt<<"]" ;
               }
               dcnt++;
            //}
            ++mtch; 
          }
          vector<uint32_t> obs_tids_vec(obs_tids.size());
          list<uint32_t>::const_iterator it = obs_tids.begin();
          const list<uint32_t>::const_iterator is = obs_tids.end();
          for(unsigned i =0; it != is; ++it, ++i) {
             obs_tids_vec[i] = *it;
          } 
          CmpDepth1 cd(dmap);
          sort(obs_tids_vec.begin(),obs_tids_vec.end(),cd);
          int last_depth = -1; 
          for(unsigned i = 0; i < obs_tids_vec.size(); ++i) {
             const uint32_t tid=obs_tids_vec[i];
             const int depth = (*dmap.find(tid)).second;
             if( depth == 0 ) {
               break;
             } 
             if( last_depth == depth || last_depth == -1 ) {
                  vector<uint32_t> path;
                  TaxTree& tax_tree_tmp = const_cast<TaxTree&>(tax_tree);
                  tax_tree_tmp.getPathToRoot(tid,path);
                  for(unsigned p = 0; p < path.size(); ++p) {
                     label_vec[pos].second.insert( make_pair(path[p],1) );
                     is_forward[pos] = !revStrnd;
                     if( tax2idx.find(tid) == tax2idx.end() ) {
                        const unsigned idx = taxid_lst.size();
                        tax2idx[tid] = idx;
                        idx2tax[idx] = tid;
                        taxid_lst.push_back(tid);
                     }
                  }
                
             } else {
               break; 
             }             
          } 
          if( verbose ) cout<<"num taxids: "<<taxid_lst.size();
          if(dcnt > 0 && verbose ) cout<<" end k-mer lookup"<<endl;
          if(mtch == 0 && verbose ) cout<<" no k-mer matches "<<endl;
          delete h;
       }
     }
   }

   void proc_line(const TaxTree& tax_tree, int ri_len, string &line, int k_size, TaxTable *table, ofstream &ofs, float threshold, const ScoreOptions& sopt, uint16_t max_count) {

     //int thread_num = omp_get_thread_num();
     assert(ri_len >= 0 && ri_len <= (signed)line.length());

       for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
          string rval;
          if( ri_len < (signed)line.length() ) {
            rval = line.substr(chk,ri_len);
          } else {
            rval = line;
          }
          const unsigned read_len = rval.length();
          if( static_cast<signed>(read_len) < ri_len ) {
            break;
          }  
          if( (signed)read_len < k_size ) {
            if( chk == 0 ) /// read too short
               ofs<<"read_label 0 -1 ReadTooShort "<<endl; 
            break;
          } 
          vector<label_info_t> label_vec(read_len-k_size+1);
          vector<bool> is_forward(read_len-k_size+1,false);
          list<uint32_t> taxid_lst; 
          hmap_t tax2idx, idx2tax;

          retrieve_kmer_labels(table, rval, k_size, label_vec, taxid_lst, tax2idx, idx2tax, false, threshold, is_forward, sopt._imap, tax_tree, max_count);

          string rev_cmp(rval.size(),'\0');
          assert(rval.size()>0);
          for (int j=(rval.size()-1); j >= 0; --j) {
            const int bidx = (rval.size()- 1) - j;
            rev_cmp[j] = getComp(rval[bidx]);
          }
          retrieve_kmer_labels(table, rev_cmp, k_size, label_vec, taxid_lst, tax2idx, idx2tax, true, threshold, is_forward, sopt._imap, tax_tree, max_count);

          //horribly slow probably
          //find LCA for retrieved nodes to make sure it's part of the set
          if( !taxid_lst.empty() ) {  
#if 0
            if( taxid_lst.size() > 1 ) {
               list<uint32_t>::const_iterator it = taxid_lst.begin();
               const list<uint32_t>::const_iterator is = taxid_lst.end();
               set<uint32_t> taxset;   

               for(; it != is; ++it) {
                  taxset.insert(*it);
               }
               const TaxNode* lca = tax_tree.lowestCommon(taxset, thread_num);
               if(!lca) {
                  if(verbose) {
                     cerr<<"Error did not get lca for: "<<line<<endl;
                     it = taxid_lst.begin();
                     for(; it != is; ++it) {
                        cerr<<"taxid="<<*it<<endl;
                     }
                  }
               } else {
                  const int lca_taxid = lca->id();
                  if( taxset.find(lca_taxid) == taxset.end() ) {
                     assert( tax2idx.find(lca_taxid) == tax2idx.end() );
                     const unsigned idx = taxid_lst.size();
                     tax2idx[lca_taxid] = idx;
                     idx2tax[idx] = lca_taxid;
                     taxid_lst.push_back(lca_taxid);
                  }
               }
            }
#endif
            construct_labels(tax_tree,label_vec,taxid_lst,tax2idx,idx2tax,ofs,k_size,sopt,is_forward); 
         } else {
            ofs<<"read_label 0 -1 ReadTooShort"<<endl; 
         }
       }
      
   }


   size_t *split_file(int n_threads, ifstream &file)
   {
     size_t * arr = new size_t[n_threads+1];

     arr[0] = 0;

     file.seekg(0, ios::end);
  size_t end = file.tellg();

  for (size_t i = 1; (signed)i<n_threads; i++)  {

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
  cout << "[-v <pct cutoff>]  -c <tax tree file> -k <kmer size> [-z:verbose] [-r:restore persistent db] [-a:ascii input format][-h:max_tid_count\n";
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


   string kmer_db_fn, query_fn, ofname, ofbase, tax_tree_fn, tax_tree_options, gid_to_tid_fn, depth_file, rand_hits_file;
   hmap_t imap;	
   ScoreOptions sopt(imap);

   char * dumpfile = NULL;


   size_t mmap_size = 0;
   
   int max_reads =0;
   uint16_t max_count = ~0;

   while ((c = getopt(argc, argv, "h:n:jb:ye:wmpk:c:v:k:i:d:l:t:s:r o:x:f:g:z:q:a ")) != -1) {
      switch(c) {

      case 'h' :
        max_count = atoi(optarg);
        if (max_count < 0) {
          max_count *= -1;
          add_root_on_kmer_drop = true;
        }
        break;
      case 's':
         mmap_size = atoi(optarg);
         mmap_size = mmap_size * (1<<30);
         cout << "Input heap size: " << mmap_size << endl;
         break;
      case 'n':
         rand_hits_file=optarg;
         break;
      case 'j':
         verbose=true;
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

     cout << "num kmers: " << taxtable->size() << endl ;

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

   loadRandHits(rand_hits_file,sopt._rand_hits);

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

     //     cout<<"Read query file: "<<query_fn<<endl;
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

       if ((pos >= arr[1+omp_get_thread_num()]) || ((signed)pos == -1)) {
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
	  proc_line(tax_tree, use_len, read_buff, k_size, taxtable, ofs, threshold,sopt, max_count);
          read_buff="";
	  read_count ++;
	  if ((signed)read_count == max_reads)
	   finished = true;

       }
     }
     ofs.close();
   }

   cout << "query time: " << clock.stop() << endl;

   return 0; 
}
