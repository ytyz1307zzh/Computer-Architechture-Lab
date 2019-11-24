#include "memory.h"

void Memory::HandleRequest(uint64_t addr, int bytes, int read,
                          char *content, int &hit, int &cycle, int &replace_cnt) {
  cycle += latency_.hit_latency;
  stats_.access_time += latency_.hit_latency;
}

