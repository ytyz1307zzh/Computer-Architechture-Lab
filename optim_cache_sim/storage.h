#ifndef CACHE_STORAGE_H_
#define CACHE_STORAGE_H_

#include <stdint.h>
#include <stdio.h>
#include <iostream>

using namespace std;

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&); \
    void operator=(const TypeName&)

struct StorageStats {
    int access_cnt;
    int miss_cnt;
    int latency;
    int replace_cnt;
    int bypass_cnt;

    StorageStats()
    {
        access_cnt = 0;
        miss_cnt = 0;
        latency = 0;
        replace_cnt = 0;
        bypass_cnt = 0;
    }
};

struct StorageLatency {
    int hit_latency;
    int bus_latency;
    StorageLatency() {}
    StorageLatency(int _hit, int _bus)
    {
        hit_latency = _hit;
        bus_latency = _bus;
    }
};

class Storage { 
    public:
	StorageStats stats_;
	StorageLatency latency_;
    Storage() {}
    ~Storage() {}

    // Sets & Gets
    void SetStats(StorageStats ss) { stats_ = ss; }
    void GetStats(StorageStats &ss) { ss = stats_; }
    void SetLatency(StorageLatency sl) { latency_ = sl; }
    void GetLatency(StorageLatency &sl) { sl = latency_; }

    void printStats(double &missrate, bool memory)
    {
        printf("access count: %d\n", stats_.access_cnt);
        printf("access time: %d\n", stats_.latency);
        if (!memory) {
            printf("miss count: %d\n", stats_.miss_cnt);
            missrate = (double)stats_.miss_cnt / (stats_.access_cnt - stats_.bypass_cnt);
            printf("miss rate: %lf\n", (double)stats_.miss_cnt / (stats_.access_cnt - stats_.bypass_cnt));
            printf("replacement count: %d\n", stats_.replace_cnt);
        }
    }

  virtual void HandleRequest(uint64_t addr, int bytes, int read,
                             char *content, int &hit, int &time, int prefetch_flag) = 0;

};

#endif //CACHE_STORAGE_H_
