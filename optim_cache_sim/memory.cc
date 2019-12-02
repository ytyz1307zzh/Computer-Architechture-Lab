#include "memory.h"

void Memory::HandleRequest(uint64_t addr, int bytes, int read,
                          char *content, int &hit, int &cycle, int prefetch_flag)
{
	// cout << "Memory prefetch flag: " << prefetch_flag << endl;
	if (!prefetch_flag)
	{
		cycle += latency_.hit_latency + latency_.bus_latency;
		stats_.latency += latency_.hit_latency + latency_.bus_latency;
		stats_.access_cnt++;
	}
}

