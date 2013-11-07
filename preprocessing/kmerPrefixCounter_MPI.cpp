#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <cassert>
#include <string>
#include <string.h>
#include <stdint.h>
#include "all_headers.hpp"
#include "Encoder.hpp"
#include <ext/hash_set>
#if HAVE_MPI == 1
#include "mpi.h"
#endif


#define HAVE_MPI 1

using namespace std;
using namespace metag;

void usage(char *argv[]) {
  cerr << "Usage:\n"
       << " -i <string>  - input fasta_fn\n"
       << " -k <int>     - kmer length\n"
       << " -o <string>  - output filename\n"
       << " -l <int>     - prefix length (as in: the string representation of the prefix)\n"
       << " -f <int> [optional] first prefix\n"
       << " -h           - print help and exit\n\n"
       << "Required arguments: i, k, o, l\n";
}

typedef uint32_t tid_t;

//void write_ascii(__gnu_cxx::hash_map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn); 
void write_ascii(map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn); 
//void write_binary(__gnu_cxx::hash_map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn);
void write_binary(map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn);

int rank = 0;
int size = 1;

int main(int argc, char *argv[]) {
  #if HAVE_MPI == 1
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  if (rank == 0) {
    cout << "sizeof(tid_t): " << sizeof(tid_t) << endl;
    cout << "invocation: ";
    for (int j=0; j<argc; j++) {
      cout << argv[j] << " ";
    }
    cout << endl;
  }  

  string input_fn, output_fn;
  int mer_len = 0;
  bool prn_help = false;
  uint64_t prefix = 0;
  int p_len = 0;
  bool ascii = false;
  int quit_early = 0;
  int first_prefix = 0;

  int count = 0;
  const string opt_string="i:k:h o:p:l:aq:f:";
  char c;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'f':
      first_prefix = atoi(optarg);
      break;
    case 'q':
      quit_early = atoi(optarg);
      break;
    case 'a':
      ascii = true;
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

  if (prn_help || count != 4) {
    if (rank == 0) {
      usage(argv);
    }
    exit(-1);
  }

  prefix = first_prefix + rank;

  uint64_t kmerid, kmerid2, rc, kmer;
  string header, seq, mer, prefix_str;
  int moveme = (mer_len*2) - (p_len*2);
  //__gnu_cxx::hash_map<uint64_t, set<tid_t> > hash;
  map<uint64_t, set<tid_t> > hash;

  ifstream in;
  if (rank == 0) {
    cout << "opening: " << input_fn << endl;
    in.open(input_fn.c_str());
    if (!in) {
      cerr << "failed to open " << input_fn << " for reading\n";
      exit(1);
    }
  }
  int s_ct = 0;
  bool finished = false;
  int sz = 0; //stop compiler complaints
  int gid = 0; //stop compiler complaints
  uint64_t bp = 0;
  StopWatch clk;
  clk.start();
  while (true) {
    //P_0 reads the headers and sequeneces
    if (rank == 0) {
      header = "";
      in>>header;
      if (header.size() == 0) {
        finished = true;
      } else {
        ++s_ct;
        if (header[0] != '>') {
          cerr << "header[0] != '>' for line number " << s_ct*2 << endl;
        }

        cerr << s_ct << " header: " << header << endl;
        assert(header[0] == '>');
        string id_s = header.substr(1);
        gid = atoi(id_s.c_str());
        in >> seq;
        cerr << "seq size: " << seq.size() << endl;
        sz = seq.size();
        bp += sz;
      }

      if (rank == 0) cerr << rank << " seq count: " << s_ct << "  bp count: " << bp << endl;

      if (finished) {
        sz = -1;
        in.close();
      } else {
        sz = seq.size();
      }

      //exit condition, for development and testing
      if (quit_early && s_ct == quit_early) {
        sz = -1;
        in.close();
      }  
    }

    //P_0 bcasts header and sequence to all others
    #if HAVE_MPI == 1
    MPI_Bcast(&sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #endif

    //exit condition
    if (sz == -1) {
      break;
    }

    #if HAVE_MPI == 1
    MPI_Bcast(&gid, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
      seq.resize(sz);
    }
    MPI_Bcast(&seq[0], sz, MPI_BYTE, 0, MPI_COMM_WORLD);
    #endif

    //hash kmers and genome IDs
    Encoder e(seq, mer_len);
//string mer1, mer2, mer3;
    while (e.next(kmerid)) {
      rc = Encoder::rc(kmerid, mer_len);
      kmer = kmerid < rc ? kmerid : rc;

/*
Utils::decode(kmerid, mer_len, mer1);
Utils::decode(rc, mer_len, mer2);
Utils::decode(kmer, mer_len, mer3);
cout << mer1 << endl;
cout << mer2 << endl;
cout << mer3 << endl << endl;
*/

      kmerid2 = kmer >> moveme;

      if (prefix == kmerid2) {
        hash[kmer].insert(gid);
      }
    }

  }

  if (rank == 0) {
    cerr << rank << " kmer count:                " << hash.size() << endl;
    cerr << "time to count kmers: " << clk.stop() << endl;
    cerr << "(stats are for P_0 only)\n";
  }  

  clk.reset();
  clk.start();
  char buf[1024];
  if (ascii) {
    sprintf(buf, "%s.%d", output_fn.c_str(), (int)prefix);
    write_ascii(hash, mer_len, buf);
  } else {
    sprintf(buf, "%s.%d", output_fn.c_str(), (int)prefix);
    if (rank == 0) cout << "writing to: " << buf << endl;
    write_binary(hash,  mer_len, buf);
  }
  if (rank == 0) {
    cerr << endl << "time to write kmers: " << clk.stop() << endl;
  }

  #if HAVE_MPI == 1
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  #endif
}


void write_ascii(map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn) {
//void write_ascii(__gnu_cxx::hash_map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn) {
    uint64_t ct = 0;
    ofstream out(fn);
    assert(out);
    string mer;
    //for (__gnu_cxx::hash_map<uint64_t, set<tid_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
    for (map<uint64_t, set<tid_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
      //out << t->first << endl;
      ++ct;
      Utils::decode(t->first, mer_len, mer);
      out << mer << " ";
      for (set<tid_t>::iterator t2 = t->second.begin(); t2 != t->second.end(); t2++) {
        out << *t2 << " ";
      }
      out << endl;
    }
    out.close();
    cerr << "num written: " << ct << endl;
}

//void write_binary(__gnu_cxx::hash_map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn) {
void write_binary(map<uint64_t, set<tid_t> > &hash, int mer_len, const char *fn) {
    KmerFileMetaData metadata;
    metadata.setSize(hash.size());
    metadata.setKmerLength(mer_len);
    metadata.setDefaultDataStart();
    if (rank == 0) cout << "opening: " << fn << endl;
    FILE *out = fopen(fn, "wb");
    assert(out);
    metadata.write(out);
    uint64_t sanity = ~0;
    tid_t  id;
    int  count;
    uint64_t ct = 0;
    uint64_t kmer;
    string mer;
    //for (__gnu_cxx::hash_map<uint64_t, set<tid_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
    for (map<uint64_t, set<tid_t> >::iterator t = hash.begin(); t != hash.end(); t++) {
      ++ct;
      kmer = t->first;
      assert(fwrite(&kmer, 8, 1, out) == 1);
      count = t->second.size();
      assert(fwrite(&count, 4, 1, out) == 1);
      for (set<tid_t>::iterator t2 = t->second.begin(); t2 != t->second.end(); t2++) {
        id = *t2;
        assert(fwrite(&id, 4, 1, out) == 1);
      }
      if (ct % 1000 == 0) {
        assert(fwrite(&sanity, 8, 1, out) == 1);
      }
    }
    fclose(out);
    if (rank == 0) cerr << "num written: " << ct << " of " << hash.size() << endl;
}
