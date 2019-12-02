#include "cache.h"
#include "def.h"
#include <cstring>
#include <iostream>
#include <cmath>

using namespace std;

extern double hit_rate1;

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

//prefetch_flag = 1: prefetch
void Cache::HandleRequest(uint64_t addr, int bytes, int read,
                          char *content, int &hit, int &cycle, int prefetch_flag) {
    // cout << "CACHE L" << level << " prefetch flag: " << prefetch_flag << endl;
    if (!prefetch_flag) 
        stats_.access_cnt++;

    int s = log2(set_num);
	int b = log2(block_num);
    uint64_t tag = getTag(addr, s, b);
    int set_idx = getSet(addr, s, b);
    int block_idx = getBlock(addr, s, b);

    int hit_idx = -1;
    for (int i = 0; i < associativity; i++) {
        if (set[set_idx].line[i].valid == 1 && set[set_idx].line[i].tag == tag)
            hit_idx = i;
    }

    if (!prefetch_flag) {
        cycle += latency_.bus_latency + latency_.hit_latency;
        stats_.latency += latency_.bus_latency + latency_.hit_latency;
    }
    // Miss?
    if (hit_idx == -1) {
        if (!prefetch_flag)
            stats_.miss_cnt++;
        ReplaceAlgorithm (addr, bytes, read, content, hit, cycle, prefetch_flag);
    }
    else {
        int s = log2(set_num);
        int b = log2(block_num);
        int setIndex = getSet(addr, s, b);
        int blockIndex = getBlock(addr, s, b);
        if (!prefetch_flag) {
            hit++;
            set[setIndex].line[hit_idx].last_access = stats_.access_cnt;
        }

        //write
        if (!read) {
            memcpy(set[setIndex].line[hit_idx].block + blockIndex, content, bytes);
            // write back
                set[setIndex].line[hit_idx].dirty = 1;
        }
        if (read)
            memcpy(content, set[setIndex].line[hit_idx].block + blockIndex, bytes);
    }
    // Prefetch
    if (!prefetch_flag && !no_optim)
        PrefetchAlgorithm(addr, hit_idx);
    return;
}

void Cache::ReplaceAlgorithm(uint64_t addr, int bytes, int read,
                             char *content, int &hit, int &cycle, int prefetch_flag)
{

    int s = log2(set_num);
    int b = log2(block_num);
    uint64_t tag = getTag(addr, s, b);
    int set_idx = getSet(addr, s, b);
    int block_idx = getBlock(addr, s, b);

    int victim = ReplaceDecision(set_idx);
    Line &victim_line = set[set_idx].line[victim];

    // the block is dirty: write to lower memory
    if (victim_line.dirty) {
        uint64_t write_addr = (victim_line.tag << (s + b)) + (set_idx << b);
        lower_ -> HandleRequest(write_addr, block_num, FALSE, victim_line.block, hit, cycle, prefetch_flag);
    }

    //request lower level to get the data
    uint64_t write_addr = addr - block_idx;
    char* newBlock = victim_line.block;
    lower_ -> HandleRequest(write_addr, block_num, TRUE, newBlock, hit, cycle, prefetch_flag);

    if (!prefetch_flag)
        stats_.replace_cnt++;

    victim_line.valid = 1;
    victim_line.reach = 1;
    victim_line.tag = tag;
    victim_line.last_access = stats_.access_cnt;

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

int Cache::ReplaceDecision(int set_idx)
{
    bool zero_flag = false;
    for (int i = 0; i < associativity; i++) {
        if (set[set_idx].line[i].reach == 0) {
            zero_flag = true;
            break;
        }
    }
    // every line has a R=1, turn to LRU
    if (!zero_flag || no_optim) {
        int victim = -1;
        int LRU_cnt = INF;
        for (int i = 0; i < associativity; i++) {
            if (set[set_idx].line[i].last_access < LRU_cnt) {
                LRU_cnt = set[set_idx].line[i].last_access;
                victim = i;
            }
            if (!no_optim)
                set[set_idx].line[i].reach = 0;
        }
        return victim;
    }
    // workset policy
    int victim = -1;
    int LRU_cnt = INF;
    for (int i = 0; i < associativity; i++) {
        if (set[set_idx].line[i].reach == 1) {
            set[set_idx].line[i].last_access = stats_.access_cnt;
            set[set_idx].line[i].reach = 0;
        }
        // R bit is zero
        else {
            // not in workset, choose this line as victim
            if (set[set_idx].line[i].last_access < stats_.access_cnt - workset_length)
                return i;
            if (set[set_idx].line[i].last_access < LRU_cnt) {
                LRU_cnt = set[set_idx].line[i].last_access;
                victim = i;
            }
        }
    }
    return victim;
}

void Cache::PrefetchAlgorithm(uint64_t addr, int hit_idx)
{
    char tmp_content[64];
    int tmp_hit = 0;
    int tmp_cycle = 0;

    for (int i = 1; i < prefetch_num; i++)
        HandleRequest(addr + 64 * i * prefetch_stride, 0, TRUE, tmp_content, tmp_hit, tmp_cycle, TRUE);
    if (hit_idx != -1 && prefetch_num > 0)
        HandleRequest(addr + 64 * prefetch_num * prefetch_stride, 0, TRUE, tmp_content, tmp_hit, tmp_cycle, TRUE);
}

//initialize the cache
void Cache::SetConfig(int _size, int _asso, int _set, int _through, int _stride,
                      int _allocate, int _prefetch, int _workset, bool _optim)
{
    size = _size;
    associativity = _asso;
    set_num = _set;
    write_through = _through;
    write_allocate = _allocate;
    block_num = size / (set_num * associativity);
    prefetch_num = _prefetch;
    prefetch_stride = _stride;
    workset_length = _workset;
    no_optim = _optim;
    set = new Set[set_num];
    for (int i = 0; i < set_num; i++) {
        set[i].line = new Line[associativity];
        memset(set[i].line, 0, associativity * sizeof(Line));
        for (int j = 0; j < associativity; j++)
            set[i].line[j].block = new char[block_num];
		set[i].replace_cnt = 0;
    }
}
