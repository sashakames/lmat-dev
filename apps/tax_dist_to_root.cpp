#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cassert>
#include <vector>
#include "all_headers.hpp"
#include <stack>
#include <queue>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::set;
using std::fstream;
using std::ofstream;
using std::ostream;

using namespace metag;

StopWatch cl;
int ctt = 0;

int main(int argc, char* argv[]) 
{
   char c = '\0';
   string taxtree_fn;
   string gid_to_tid_fn;
   bool prn_help=false;
   const string opt_string="t:h:g:";   
   while ((c = getopt(argc, argv, opt_string.c_str())) != -1) {
      switch(c) {
      case 'g':
      cout << optarg << endl;
        gid_to_tid_fn = optarg;
        break;
      case 't':
         taxtree_fn = optarg;
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
      exit(0);
   }

   cout << "Start taxtree load\n";

   //dah I'm hardcoding to show what files to use; should be
   //    changed to use cmd line parameters
   TaxTree tax_tree(taxtree_fn.c_str(), gid_to_tid_fn.c_str(), 1);
   //dah: no longer needed; kpath_taxonomy.dat contains both parent and children for each node
   //tax_tree.findChildren();
   //

   TaxTree::const_iterator it = tax_tree.begin();
   TaxTree::const_iterator is = tax_tree.end();
   for(; it != is; ++it) {
      vector<uint32_t> path;
      const uint32_t taxid = (*it).first;
      tax_tree.getPathToRoot(taxid,path);
      cout<<taxid<<" "<<path.size()<<endl;
   }
   return 0;
   
}
