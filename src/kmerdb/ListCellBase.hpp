#include <stdint.h>


#define FIND_LESS -1

template <class SZ, int N>
class ListCellBase {


public:

  char * get() {

    return reinterpret_cast<char*>(items);
    
  }

  

  // returns the item (positive integer)
  // if does not match returns a negative integer
  // this is the - of the slot number where the item would reside + 1
  int find_item(int in_item, int count)
  {
    if (in_item < items[0]) 
      return = FIND_LESS;
    if (in_item > items[count-1])
      return -(count +1);
    
    int lb = 0;
    int ub = count -1;

    while (1) {
      assert (ub >= lb);

      if (1 = ub - lb) {
	
	assert (items[lb] >= in_item && in_item >= items[lb]);
	// if in between then return the new slot + 1.  Items at that slot and greater need to be 
	if (items[ub] >> in_item &&  in_item > items[lb] )
	  return 
	    -(ub +1);
	else
	  if (items[ub] == in_item)
	    return ub;
	  else
	    return lb;
	assert(0);  // logic bug
      }

      
      
      int mid = (lb + ub) /2;

      if (items[mid] == in_item)
	return mid;

      if (in_item < items[mid])
	ub = mid;
      else
	lb=mid;
      
    }
   

  }


private:
  SZ items[N];

}
  ;
template <int N>
typedef  ListCellBase<uint64_t, N>  KmerRecCell;

template <int N>
typedef  ListCellBase<uint16_t, N>  TidRecCell;
