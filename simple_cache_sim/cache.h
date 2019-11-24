#ifndef CACHE_CACHE_H_
#define CACHE_CACHE_H_

#include <stdint.h>
#include "storage.h"

struct Line{
    int valid;
    uint64_t tag;
    int last_access;
    int dirty;
    char *block;
};

struct Set{
    Line *line;
};

class Cache: public Storage {
    public:
    Cache() {}
    ~Cache() {}

    // Sets & Gets
    void SetConfig(int _size, int _asso, int _set, int _through, int _allocate, int _debug);
    void SetLower(Storage *ll) { lower_ = ll; }
    // Main access process
    void HandleRequest(uint64_t addr, int bytes, int read,
                         char *content, int &hit, int &time, int &replace_cnt);


    private:
    int size;
    int associativity;
    int set_num;
    int block_num;
    int write_through;
    int write_allocate;
    int debug;
    int access_count;

    int ReplaceDecision(uint64_t addr, int &index);
    void ReplaceAlgorithm(uint64_t addr, int bytes, int read,
                          char *content, int &hit, int &cycle, int &replace_cnt);

    Storage *lower_;
    Set *set;
    DISALLOW_COPY_AND_ASSIGN(Cache);
};

#endif //CACHE_CACHE_H_
