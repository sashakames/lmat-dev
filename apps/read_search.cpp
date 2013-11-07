#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "kencode.hpp"
#include "all_headers.hpp"

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace kencode_ns;
using namespace metag;

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

void print_mer_counts(const KmerDB& h, const string& read, const size_t mer_len, vector<bool>& pa_vec, bool revStrnd)
{
//dah debug:
  kencode_c ken(mer_len);
  size_t read_len = read.length();
  uint64_t kmer_id; 
  bool bad_last_kmer=true;
  unsigned next_valid_pos = 0;
  for(size_t j = 0; j < read.size(); ++j) {
    if(!isBase(read[j])) {
        bad_last_kmer=true;
        next_valid_pos = (j + mer_len);
        continue; 
    }
    if( bad_last_kmer && (j >= mer_len - 1) && j >= next_valid_pos ) { 
        const string nkmer=read.substr((j-mer_len)+1,mer_len);
        kmer_id = ken.kencode(nkmer);
        bad_last_kmer=false;
    } else if( !bad_last_kmer) {
        kmer_id = ken.kencode(read[j]);
        bad_last_kmer=false;
    }
    if(!bad_last_kmer) {
       KmerDB& tfind = const_cast<KmerDB&>(h);
       if (tfind.lookup(kmer_id)) {
         size_t pos = j - mer_len + 1;
         if (revStrnd) {
            pos = read_len - j - 1;
         }
         pa_vec[pos] = true;
       }
    }
  }
}

int main(int argc, char* argv[]) 
{

   char c = '\0';
   int k_size=-1, ri_len = -1;
   string kmer_db_fn, query_fn;
   while ((c = getopt(argc, argv, "k:i:d:l:")) != -1) {
      switch(c) {
      case 'l':
         ri_len = atoi(optarg);
         break;
      case 'i':
         query_fn = optarg;
         break;
      case 'd':
         kmer_db_fn = optarg;
         break;
      default:
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }
   cout << "Start kmer DB load\n";
   KmerDB table;
   table.ingest(kmer_db_fn.c_str());
   cout << "End kmer DB load\n";

   k_size = table.get_kmer_length();
   kencode_c ken(k_size);

   ifstream ifs(query_fn.c_str());
   string line;
   /* ASSUME READ FITS ON ONE LINE FOR NOW!! */
   while (getline(ifs, line))   {
    if (line[0] == '>') {
      continue;
    }
    if( ri_len == -1) {
      // if not specified use first read as the fixed length
      ri_len = line.length();
    } 
    cout<<"read: "<<line<<endl;
    //split read into multiples of "ri_len" (read interval length)
    for(unsigned chk = 0; (chk+ri_len) <= line.length(); chk += (ri_len+1)) {
       string rval = line.substr(chk,ri_len);
       const unsigned read_len = rval.length();
       if( static_cast<signed>(read_len) < ri_len ) {
         break;
       }  
       vector<bool> pa_vec(read_len - k_size + 1);

       print_mer_counts(table, rval, k_size, pa_vec, false);

       string rev_cmp(rval.size(),'\0');
       assert(rval.size()>0);
       for (int j=(rval.size()-1); j >= 0; --j) {
         const int bidx = (rval.size()- 1) - j;
         rev_cmp[j] = getComp(rval[bidx]);
       }
       print_mer_counts(table, rev_cmp, k_size, pa_vec, true);
       unsigned kmer_cnt = 0;
       for(unsigned i = 0; i < read_len; ++i) {
          if(pa_vec[i]) {
            ++kmer_cnt;
          }
       }
       const float fval = static_cast<float>(kmer_cnt) / static_cast<float>(read_len);
       cout<<fval<<" "<<kmer_cnt<<" "<<read_len<< endl;
     }
   }
   return 0; 
}
