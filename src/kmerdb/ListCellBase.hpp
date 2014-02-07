#include <stdint.h>

#define FIND_LESS -1

template <class SZ, int N>
class ListCellBase {


public:

  char * get() {

    return reinterpret_cast<char*>(m_items);
    
  }


    bool insert_item(SZ in_item, int count) {
        
        
        char *buf[sizeof(SZ) * 2048];
        
        int pos = find_item(in_item, count);
        
        if (pos >= 0)
            return false;
        pos = (-1 - pos);
        
        if (pos < count)
        {
            char *addr = get() + (sizeof(SZ) * (pos+1));
            memcpy(buf, addr, (sizeof(SZ) * (count-pos)));
            memcpy(get()+(sizeof(SZ) * pos+1), buf, (sizeof(SZ) * (count-pos)) );
                   
        }
        
        m_items[pos] = in_item;
        
    }
    

  // returns the item (positive integer)
  // if does not match returns a negative integer
  // this is the - of the slot number where the item would reside + 1
  int find_item(SZ in_item, int count)
  {
    if (in_item < m_items[0]) 
      return FIND_LESS;
    if (in_item > m_items[count-1])
      return -(count +1);
    
    int lb = 0;
    int ub = count -1;

    while (1) {
      assert (ub >= lb);

      if (1 == ub - lb) {
	
	assert (m_items[lb] >= in_item && in_item >= m_items[lb]);
	// if in between then return the new slot + 1.  M_Items at that slot and greater need to be 
	if (m_items[ub] >> in_item &&  in_item > m_items[lb] )
	  return 
	    -(ub +1);
	else
	  if (m_items[ub] == in_item)
	    return ub;
	  else
	    return lb;
	assert(0);  // logic bug
      }

      
      
      int mid = (lb + ub) /2;

      if (m_items[mid] == in_item)
	return mid;

      if (in_item < m_items[mid])
	ub = mid;
      else
	lb=mid;
      
    }
   

  }


  SZ m_items[N];

};



