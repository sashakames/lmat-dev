#include <unistd.h>
#include <cassert>
#include <cstdio>


// return total program size and resident set size
void memory(int& size, int& rss) {
  static int ps = getpagesize();
  FILE* file = fopen("/proc/self/statm", "r");
  assert(file);
  int count = fscanf(file, "%d%d", &size, &rss);
  assert(count == 2);
  size *= ps;
  rss *= ps;
  fclose(file);
}
