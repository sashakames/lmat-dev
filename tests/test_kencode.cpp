#include <unistd.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdint.h>
#include "kencode.hpp"

using namespace std;

using namespace kencode_ns;

int main(int argc, char *argv[]) 
{

  int KSize=6;                  // "k" - length of the k-mer
  
  string str, str1;
  uint64_t kmerid;
  int pos=0;
  int c;


  while ((c = getopt(argc, argv, "s:")) != -1) {
    switch(c) {
    case 's':
      KSize = atoi(optarg);
      break;
    default:
      cout << "kmer [kmer-size], default is " << KSize << "\n";
    }
  }

  kencode_c ken(KSize);
  str1="                                                                          ";
  while (getline(cin, str)) {
    if (str.empty()) continue; 
    if (str.size() < (unsigned) KSize) {
	cout << "String must be at least " << KSize << endl;
	continue;
      }
    
    
    kmerid = ken.kencode(str);
    printf("string, id %s %llx\n", str.c_str(), kmerid);
    ken.kdecode(kmerid, str1);
    cout << "'" << str1 << "'" <<endl;
    for (pos=KSize; pos < (signed) str.size(); pos +=1) {
      kmerid = ken.kencode(str[pos]);
      printf("string, id %c %llx\n", str[pos], kmerid);
      ken.kdecode(kmerid, str1);
      cout << str1 << endl;
    }

  }
 
  return(0);
}
