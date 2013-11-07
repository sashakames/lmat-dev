#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include <cstdlib>
#include <map>
#include <algorithm>

using namespace std;

void
usage();

struct cmp :
  public binary_function<const pair<int,int>,const pair<int,int>,bool >  {
    bool operator()(const pair<int,int> &a, const pair<int,int>&b) {
      return a.first > b.first;
    }
};

int main(int argc, char* argv[]) {

  cout << "invocation: ";
  for (int j=0; j<argc; j++) {
    cout << argv[j] << " ";
  }
  cout << endl;

  char c = '\0';
  bool prn_help = false;
  string output_fn, histo_fn, read_fn;
  size_t mer_len = 0;
  const string opt_string="o:d:k:g:h q:n:";
  int count = 0;
  int quit_early = 0;
  int n = 0;
  while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
    switch (c) {
    case 'n':
      n = atoi(optarg);
      ++count;
      break;
    case 'q':
      quit_early = atoi(optarg);
      break;
    case 'o':
      output_fn = optarg;
      break;
    case 'd':
      read_fn = optarg;
      ++count;
      break;
    case 'g':
      histo_fn = optarg;
      break;
    case 'k':
      mer_len = atoi(optarg);
      ++count;
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
  if ( prn_help || count != 3 ) {
    usage();
    exit(0);
  }

  ofstream saveme;
  if (output_fn.size()) {
    saveme.open(output_fn.c_str());
    assert(saveme);
  }

  ifstream in(read_fn.c_str());
  assert(in);
  string header, seq;
  int ct = 0;
  int total = 0;
  int good = 0;
  map<int, int> mp;
  while (getline(in, header) > 0) {
    ++ct;
    if (ct == quit_early) break;
    assert(getline(in, seq) > 0);
    assert(header[0] == '>');
    assert(seq[0] != '>');
    for (size_t j=0; j<seq.size(); j++) {
      char c = tolower(seq[j]);
      if (! (c == 'a' || c == 'c' || c =='g' || c == 't' || c == 'n')) {
        cerr << "unknown character in read " << ct << " not a,c,g,t or n\n";
        cout << seq << endl;
        exit(1);
      }
    }
    size_t idx = -1;
    size_t last = 0;
    int num = 0;
    while (true) {
      ++idx;
      last = idx;
      idx = seq.find('N', idx);
      if (idx == string::npos) {
        if (seq.size() - last >= mer_len) {
          num += seq.size() - last - mer_len + 1;
        }
        break;
      }
      if (idx - last >= mer_len) {
          num += idx - last - mer_len + 1;
        }
    }
    if (num>= n) {
      ++good;
      if (saveme) {
        saveme << header << endl << seq << endl;
      }
    }
    if (mp.find(num) == mp.end()) {
      mp[num] = 0;
    }
      mp[num] += 1;
    ++total;

  }
  cout << "read count: " << total << endl;
  cout << "would keep:  " << good << " (" << 100.0*good/total<<"%)"<<endl;

  if (saveme) {
    saveme.close();
  }

  //optionally write histogram data
  if (histo_fn.size()) {
    ofstream out(histo_fn.c_str());
    assert(out);
    vector<pair<int,int> > v;
    for (map<int,int>::const_iterator t = mp.begin(); t != mp.end(); t++) {
      v.push_back(pair<int,int>(t->first, t->second));
    }
    sort(v.begin(), v.end(), cmp());
    for (size_t j=0; j<v.size(); j++) {
      out << v[j].first<<" "<<v[j].second<<endl;
    }
    out.close();
  }


  in.close();
  return 0;
}


void
usage(){
  cout<<"Usage: \n"
      "  -d <string>   - input fasta filename                                    [required]\n"
      "  -k <int>      - kmer length                                             [required]\n"
      "  -n <int>      - min number of kmers for considering read as usable      [required]\n"
      "  -o <string>   - output fn;                                              [optional]\n"
      "  -g <string>   - write histogram: number of usable kmers vs. read count  [optional]\n"
      "  -q <int>      - quit after 'q' sequences; for debugging and testing     [optional]\n"
      "  -h            - print help and exit                                     [optional]\n" 
      "\n"
      "Assumes: each read sequence is on a single line\n";
}
