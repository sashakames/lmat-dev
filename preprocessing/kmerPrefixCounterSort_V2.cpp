#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <cassert>
#include <string>
#include <cstdlib>
#include <string.h>
#include <stdint.h>
#include "all_headers.hpp"
#include "Encoder.hpp"
#include <ext/hash_set>


using namespace std;
using namespace metag;


void usage(char *argv[]) {
  cerr << "Usage:\n"
       << " -i <string>  - input fasta_fn\n"
       << " -k <int>     - kmer length\n"
       << " -o <string>  - output filename\n"
       << " -l <int>     - prefix length (as in: the string representation of the prefix)\n"
       << " -m <int>     - max memory, in G, available for vector allocation;\n"
       << "                should be a few G less than max available; default\n"
       << "                is 24G, which is sufficient for 1.5G kmers\n"
       << " -f <int>     - prefix\n"
       << " -a           - output kmers in ascii format (useful for testing)\n"
       << " -h           - print help and exit\n\n";
}

std::ifstream::pos_type filesize(const char* filename) {
  std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
  assert(in);
  in.seekg(0, std::ifstream::end);
  std::ifstream::pos_type r = in.tellg();
  in.close();
  return r;
}


struct P {
  P(uint64_t kmer, uint32_t tid) : m_kmer(kmer), m_tid(tid) {}
  P() {}
  uint64_t m_kmer;
  uint32_t m_tid;
};





#define ENCODE(t, c, k) \
switch (c) { \
case 'a': case 'A': t = 0; break; \
case 'c': case 'C': t = 1; break; \
case 'g': case 'G': t = 2; break; \
case 't': case 'T': t = 3; break; \
default: k = 0; continue; \
}

/*
 *    Given a sequence of nucleotide characters,
 *    break it into canonical k-mers in one pass.
 *    Nucleotides are encoded with two bits in
 *    the k-mer. Any k-mers with ambiguous characters
 *    are skipped.
 *    str:  DNA sequence (read)
 *    slen: DNA sequence length in nucleotides
 *    klen: k-mer length in nucleotides
 **/
int moveme;
uint64_t prefix = 0;

static
void seq_lookup(const char *str, int slen, int klen, int gid, vector<P> &data, uint64_t &data_idx) {
  int j; /* position of last nucleotide in sequence */
  int k = 0; /* count of contiguous valid characters */
  int highbits = (klen-1)*2; /* shift needed to reach highest k-mer bits */
  uint64_t mask = ((uint64_t)1 << klen*2)-1; /* bits covering encoded k-mer */
  uint64_t forward = 0; /* forward k-mer */
  uint64_t reverse = 0; /* reverse k-mer */

  for (j = 0; j < slen; j++) {
    register int t;
    ENCODE(t, str[j], k);
    forward = ((forward << 2) | t) & mask;
    reverse = ((uint64_t)(t^3) << highbits) | (reverse >> 2);
    if (++k >= klen) {
      /* kmer_lookup(kmer); do k-mer lookup here... */
        uint64_t kmer = forward < reverse ? forward : reverse;
        uint64_t kmer2 = kmer >> moveme;
        if (prefix == kmer2) {
        //P_global.m_kmer = kmer;
        //P_global.m_tid = gid;
        //data.push_back(P_global);
        //data.push_back(P(kmer, gid));
        //cout << ">>>>>>>> " << data_idx << endl;
        data[data_idx].m_kmer = kmer;
        data[data_idx].m_tid = gid;
        ++data_idx;
        }
      

      /* zero based position of forward k-mer is (j-klen+1) */
      /* zero based position of reverse k-mer is (slen-j-1) */
    }
  }
}

struct  myclass {
  bool operator() (P a, P b) {
    if (a.m_kmer != b.m_kmer) {
      return a.m_kmer < b.m_kmer;
    }
    return a.m_tid < b.m_tid;
  }
} mysort;



uint64_t write_binary(const vector<P> &v, int mer_len, const char *fn);
uint64_t write_ascii(const vector<P> &v, int mer_len, const char *fn);

int main(int argc, char *argv[]) {
  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  uint64_t data_idx = 0;

  string input_fn, output_fn;
  int mer_len = 0;
  bool prn_help = false;
  int p_len = 0;
  int quit_early = 0;
  bool ascii = false;
  uint64_t mem = 24;

  int count = 0;
  const string opt_string="a i:k:h o:p:l:q:f:m:";
  char c;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'm' :
      mem = atoi(optarg);
      break;
    case 'a' :
      ascii = true;
      break;
    case 'f':
      ++count;
      prefix = atoi(optarg);
      break;
    case 'q':
      quit_early = atoi(optarg);
      break;
    case 'i':
      input_fn = optarg;
      ++count;
      break;
    case 'k':
      mer_len = atoi(optarg);
      ++count;
      break;
    case 'h':
      prn_help = false;
      break;
    case 'o':
      ++count;
      output_fn = optarg;
      break;
      break;
    case 'l':
      ++count;
      p_len = atoi(optarg);
      break;
    default:
      cerr << "Unrecognized option: "<<c<<", ignore."<<endl;
      prn_help = true;
      break;
    }
  }

  if (prn_help || count != 5) {
    usage(argv);
    exit(-1);
  }


  string header, seq, mer, prefix_str;
  moveme = (mer_len*2) - (p_len*2);
  //vector<P> v;
  mem *= 1000000000;
  mem /= sizeof(P);
  cerr << "allocating: " << mem << " slots= " << mem*sizeof(P) << " bytes\n";
  vector<P> data;
//  data.reserve(mem);
  data.resize(mem);

  ifstream in;
  cout << "opening: " << input_fn << endl;
  in.open(input_fn.c_str());
  if (!in) {
    cerr << "failed to open " << input_fn << " for reading\n";
    exit(1);
  }

  int s_ct = 0;
  int gid = 0; //stop compiler complaints
  StopWatch clk;
  clk.start();
  StopWatch c2;
  while (true) {
    c2.start();
    header = "";
    in>>header;
    if (header.size() == 0) {
      break;
      in.close();
    } else {
      ++s_ct;

      //sanity check: each sequence must bo contained in a single line
      if (header[0] != '>') {
        cerr << "header[0] != '>' for line number " << s_ct*2 << endl;
      }
      assert(header[0] == '>');
      string id_s = header.substr(1);
      gid = atoi(id_s.c_str());
      in >> seq;

    //exit condition, for development and testing
    if (quit_early && s_ct == quit_early) {
      in.close();
      break;
    }

    if (data_idx + seq.size() > mem) {
      mem += seq.size();
      mem *= 2;
      cout << "reallocating! old size: " << data.size() << " new size: " << mem << endl;
      data.resize(mem);
    }
    seq_lookup(seq.c_str(), seq.size(), mer_len, gid, data, data_idx);
  }
  cout << "seq #: " << s_ct << " header: " << header << " time: " << c2.stop() << " len: " << seq.size() << endl;
  c2.reset();

  }

  cout << "allocated memory: " <<data.size() << endl;
  cout << "used memory: " <<data_idx << endl;
  data.resize(data_idx);

  double compute_time = clk.stop();
  cout << "time to encode kmers: " << compute_time << endl;
  clk.reset(); clk.start();
  sort(data.begin(), data.end(), mysort);
  double sort_time = clk.stop();
  cout << "time to sort: " << sort_time << endl;
  clk.reset(); clk.start();

  char buf[1024];
  sprintf(buf, "%s.%d", output_fn.c_str(), (int)prefix);
  uint64_t kmer_ct = 0;
  if (ascii) {
    kmer_ct = write_ascii(data,  mer_len, buf);
  } else {  
    kmer_ct = write_binary(data,  mer_len, buf);
  }  
  double write_time = clk.stop();
  cout << "time to write output: " << write_time << endl;
  cout << "\ninvocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  /*
  cout << "proc/[pid]/status:\n";
  pid_t p = getpid();
  sprintf(buf, "cat /proc/%d/status", p);
  system(buf);
  */

  cout << "kmer_ct: " << kmer_ct << " compute time: " << compute_time << " sort time: " << sort_time << " write_time: " << write_time << endl;
  cout << "SUCCESS!\n";
}


uint64_t write_binary(const vector<P> &v, int mer_len, const char *fn) {
  KmerFileMetaData metadata;
  //metadata.setSize(hash.size());
  metadata.setKmerLength(mer_len);
  metadata.setDefaultDataStart();
  cout << "writing to: " << fn << endl;
  FILE *out = fopen(fn, "wb");
  assert(out);
  metadata.write(out);
  uint64_t sanity = ~0;

  int64_t last_kmer = v[0].m_kmer;
  uint32_t last_tid = v[0].m_tid;
  uint64_t kmer_ct = 0;
  uint32_t tid_ct;
  vector<uint32_t> work;
  work.reserve(100000);

  size_t j = 0;
  while (j < v.size()) {
    ++kmer_ct;

    last_kmer = v[j].m_kmer;
    if (kmer_ct < 5) cout << "writing kmer " << last_kmer << " at offset " << ftell(out) << endl;
    assert(fwrite(&last_kmer, 8, 1, out) == 1);
    last_tid = v[j].m_tid;
    work.clear();
    work.push_back(last_tid);
    uint64_t k = j+1;
    while (true) {
      if ((int64_t)v[k].m_kmer != (int64_t)last_kmer || (int64_t)k == (int64_t) v.size()) {
        j = k;
        break;
      }
      if (v[k].m_tid != last_tid) {
        last_tid = v[k].m_tid;
        work.push_back(last_tid);
      }
      ++k;
    }
    tid_ct = work.size();
    assert(fwrite(&tid_ct, 4, 1, out) == 1);
    for (size_t k=0; k<work.size(); k++) {
      assert(fwrite(&work[k], 4, 1, out) == 1);
    }
    if (kmer_ct % 1000 == 0) {
      assert(fwrite(&sanity, 8, 1, out) == 1);
    }
  }

  metadata.setSize(kmer_ct);
  fseek(out, 0, SEEK_SET);
  metadata.write(out);

  fclose(out);
  return kmer_ct;
}

uint64_t write_ascii(const vector<P> &v, int mer_len, const char *fn) {
  cout << "writing to: " << fn << endl;
  ofstream out(fn);
  assert(out);

  string mer;

  int64_t last_kmer = v[0].m_kmer;
  uint32_t last_tid = v[0].m_tid;
  uint64_t kmer_ct = 0;
  uint32_t tid_ct;
  vector<uint32_t> work;
  work.reserve(100000);

  size_t j = 0;
  while (j < v.size()) {
    ++kmer_ct;
    last_kmer = v[j].m_kmer;
    Encoder::decode(last_kmer, mer_len, mer);
    std::transform(mer.begin(), mer.end(), mer.begin(), ::tolower);
    out << mer;
    last_tid = v[j].m_tid;
    work.clear();
    work.push_back(last_tid);
    uint64_t k = j+1;
    while (true) {
      if ((int64_t)v[k].m_kmer != (int64_t)last_kmer || (int64_t)k == (int64_t)v.size()) {
        j = k;
        break;
      }
      if (v[k].m_tid != last_tid) {
        last_tid = v[k].m_tid;
        work.push_back(last_tid);
      }
      ++k;
    }
    tid_ct = work.size();
    for (size_t k=0; k<work.size(); k++) {
      out << " " << work[k];
    }
    out << endl;
  }

  out.close();
  return kmer_ct;
}

