#include "cache.h"
#include "def.h"
#include <cstring>
#include <iostream>
#include <cmath>

using namespace std;

int log2(uint64_t x)
{
    return log(x) / log(2);
}

int getSet(uint64_t addr, int s, int b)
{
    return (addr >> b) & ((1 << s) - 1);
}

uint64_t getTag(uint64_t addr, int s, int b)
{
    return (addr >> (s + b));
}

int getBlock(uint64_t addr, int s, int b)
{
    return addr & ((1 << b) - 1);
}

void Cache::HandleRequest(uint64_t addr, int bytes, int read,
                          char *content, int &hit, int &cycle, int &replace_cnt) 
{
    access_count++;
    int s = log2(set_num);
    int b = log2(block_num);
    uint64_t tag = getTag(addr, s, b);
    int set_idx = getSet(addr, s, b);
    int block_idx = getBlock(addr, s, b);
    if (debug) {
        printf("query tag: 0x%llx\n", tag);
        printf("Memory Status:\n");
        printf("Set index: %x\n", set_idx);
    }

    int hit_idx = -1;
    for (int i = 0; i < associativity; i++) {
        if (debug)
            printf("valid: 0x%x, tag: 0x%llx\n", set[set_idx].line[i].valid, set[set_idx].line[i].tag);
        if (set[set_idx].line[i].valid == 1 && set[set_idx].line[i].tag == tag)
            hit_idx = i;
    }

    cycle += latency_.hit_latency;
    stats_.access_time += latency_.hit_latency;
    if (hit_idx == -1) {
        if (debug)
            cout << "MISS" << endl;
        replace_cnt++;
        ReplaceAlgorithm (addr, bytes, read, content, hit, cycle, replace_cnt);
    }
    else {
        // return hit & time
        hit++;
        if (debug)
            cout << "HIT" << endl;
        set[set_idx].line[hit_idx].last_access = access_count;

        // write
        if (!read) {
            memcpy(set[set_idx].line[hit_idx].block + block_idx, content, bytes);
            // write back
            if (!write_through)
                set[set_idx].line[hit_idx].dirty = 1;
            // write through
            else
                lower_->HandleRequest(addr, bytes, read, content, hit, cycle, replace_cnt);
        }
        // read
        else
            memcpy(content, set[set_idx].line[hit_idx].block + block_idx, bytes);
    }
}

void Cache::ReplaceAlgorithm(uint64_t addr, int bytes, int read,
                             char *content, int &hit, int &cycle, int &replace_cnt)
{
    // no-write allocate
    if (!read && !write_allocate) {
        lower_-> HandleRequest(addr, bytes, read, content, hit, cycle, replace_cnt);
        return;
    }

    int s = log2(set_num);
    int b = log2(block_num);
    uint64_t tag = getTag(addr, s, b);
    int set_idx = getSet(addr, s, b);
    int block_idx = getBlock(addr, s, b);

    int victim = -1;
    int LRU_cnt = INF;
    for (int i = 0; i < associativity; i++) {
        if (set[set_idx].line[i].last_access < LRU_cnt) {
            LRU_cnt = set[set_idx].line[i].last_access;
            victim = i;
        }
    }
    Line &victim_line = set[set_idx].line[victim];

    // the block is dirty: write to lower memory
    if (victim_line.dirty) {
        uint64_t write_addr = (victim_line.tag << (s + b)) + (set_idx << b);
        lower_ -> HandleRequest(write_addr, block_num, FALSE, victim_line.block, hit, cycle, replace_cnt);
    }

    //request lower level to get the data
    uint64_t write_addr = addr - block_idx;
    char* newBlock = victim_line.block;
    lower_ -> HandleRequest(write_addr, block_num, TRUE, newBlock, hit, cycle, replace_cnt);

    victim_line.valid = 1;
    victim_line.tag = tag;
    victim_line.last_access = access_count;

    // read
    if (read) {
        victim_line.dirty = 0;
        memcpy(content, victim_line.block + block_idx, bytes);
    }
    // write allocate
    else {
        victim_line.dirty = 1;
        memcpy(victim_line.block + block_idx, content, bytes);
    }
}

//initialize the cache
void Cache::SetConfig(int _size, int _asso, int _set, int _through, int _allocate, int _debug)
{
    size = _size;
    associativity = _asso;
    set_num = _set;
    write_through = _through;
    write_allocate = _allocate;
    block_num = size / (set_num * associativity);
    access_count = 0;
    debug = _debug;
    set = new Set[set_num];

    for (int i = 0; i < set_num; i++) {
        set[i].line = new Line[associativity];
        memset(set[i].line, 0, associativity * sizeof(Line));
        for (int j = 0; j < associativity; j++) {
            set[i].line[j].block = new char[block_num];
        }
    }
}
