#include <fstream>
#include <iostream>

#include <vector>
#include <all_headers.hpp>

using std::fstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;


using namespace metag;

size_t perm_bytes_allocd;

int main(int argc, char* argv[]) 
{

  string jellylist_fn, kmer_id_fn, gt_fn;

  char c;

   while ((c = getopt(argc, argv, "i:k:g:")) != -1) {
      switch(c) {
      case 'i':
	jellylist_fn = optarg;
	break;
      case 'k':
	kmer_id_fn = optarg;
	break;
      case 'g':
	gt_fn = optarg;
	break;
      default:
	cout << argv[0] << " searches jellylist files (-i) for kmer ids (-k) listed as 64-bit decimal integer.  Uses (-g) genomid_to_taxid data file." << endl;
         cout << "Unrecognized option: "<<c<<", ignore."<<endl;
      }
   }

  GenomeIdToTaxId lookup(gt_fn.c_str());


  ifstream ifs(kmer_id_fn.c_str());

  string line;

  vector<uint64_t> kid_vec;

  
  while(ifs >> line) {

    uint64_t kid;

    sscanf(line.c_str(),"%lld", &kid);

    kid_vec.push_back(kid);
  }
  cout << "loaded kmer ids: " << kid_vec.size() << endl;

  StopWatch clock;

  cout << "instantiating KmerNodes from file\n";


  KmerFileMetaData metadata;
  FILE *fp = Utils::openReadFile(jellylist_fn.c_str()); 
  metadata.read(fp);

  uint64_t kmer_count = metadata.size();
  metadata.write(); 

  int found_count = 0;

  
  clock.reset(); clock.start();
  KmerNode *w;
  uint64_t test, sanity = ~0;
    
  
  for (uint64_t j=0; j<metadata.size(); j++) {
    
    bool found = false;
    
    uint64_t kmer_id;

    assert(fread(&kmer_id, sizeof(kmer_t), 1, fp) == 1);
    uint32_t g_ct, g_id;
    assert(fread(&g_ct, sizeof(uint32_t), 1, fp) == 1);
      
    for (vector<uint64_t>::const_iterator it = kid_vec.begin(); it != kid_vec.end(); it++) {
      if (kmer_id == *it) {
	found_count++;
	found = true;
	cout << kmer_id << " ";
	  
	break;

      }
	

    }
      
    
    for (uint32_t h=0; h<g_ct; h++) {

      assert(fread(&g_id, sizeof(uint32_t), 1, fp) == 1);
	
      if (found) {


	  cout<<" "<<g_id<<":";
	  GenomeIdToTaxId::const_iterator s_iter = lookup.find(g_id);
	  if (s_iter != lookup.end()) {
	    cout<<s_iter->second;
	    
	  } 

      }


    }	  

    if (found) {
      cout<< endl;
    }
    if ((j+1) % 1000 == 0) {
      assert(fread(&test, 8, 1, fp) == 1);
      assert(test == sanity);
    }

    if (found_count == kid_vec.size())
      break;
  }
  fclose(fp);

  cout << "time to read kmers: " << clock.stop() << endl;
    
}
