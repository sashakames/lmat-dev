#ifndef _TID_CHECKS_HPP
#define _TID_CHECKS_HPP

#include <stdint.h>

#define TID_T uint32_t

/* attemp to centralize any references to harcoded tax ids here */

#define PHIX_TID 374840
//#define PHIX_TID 1217068
#define HUMAN_TID 9606

#define isPhiX(tid) (((tid==PHIX_TID)||(tid==10847)||(tid==374840)||tid==32630))? 1 : 0)

inline bool isHuman(TID_T taxid)  {
   bool res = false;
   switch (taxid) {
#if 0
   case 10000348:
   case 10000349:
   case 10000350:
   case 10000351:
   case 10000352:
   case 10000353:
   case 10000354:
   case 10000355:
   case 10000356:
   case 10000357:
   case 10000358:
   case 10000359:
   case 10000360:
   case 10000361:
   case 10000362:
   case 10000363:
   case 10000364:
   case 10000365:
   case 10000366:
   case 10000367:
   case 10000368:
   case 10000369:
   case 10000370:
   case 10000371:
   case 10000372:
   case 10000373:
#endif
   case 9606:
   case 63221: //neanderthal:
   case 741158: //denisovan
      res = true;
      break;
   default:
      res=false;
      break;
   }
   return res;
}

inline bool isEukId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 2759:
      res = true;
      break;
   default:
      res = false;
      break; 
   }
   return res;
}

inline bool isVirId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 10239:
      res = true;
      break;
   default:
      res = false;
      break;
   }
   return res;
}

inline bool isProkId(TID_T taxid) {
   bool res = false;
   switch (taxid) {
   case 2157: //Archaea
   case 2:
      res = true;
      break;
   default:
      res = false;
      break;
   }
   return res;
}


#endif // _TID_CHECKS_HPP
