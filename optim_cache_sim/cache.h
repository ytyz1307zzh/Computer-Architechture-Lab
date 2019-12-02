#ifndef CACHE_CACHE_H_
#define CACHE_CACHE_H_

#include <stdint.h>
#include "storage.h"
using namespace std;


//set结构体嵌套line结构体
struct Line
{
	int valid;
	uint64_t tag;
	int reach;
	int last_access;
	int dirty;
	char *block;
	int prefetch_flag;
};

struct Set
{
	Line *line;
	int replace_cnt;
};


class Cache : public Storage {
public:
	Cache() {}
	Cache(int _level) { level = _level; }
	~Cache() { delete set; }

	// Sets & Gets
	void SetConfig(int _size, int _asso, int _set, int _through, int _stride,
				   int _allocate, int _prefetch, int _workset, bool _optim);
	void SetLower(Storage *ll) { lower_ = ll; }

	void HandleRequest(uint64_t addr, int bytes, int read,
					   char *content, int &hit, int &time, int prefetch_flag);

private:
	int size;
    int associativity;
    int set_num;
    int block_num;
    int write_through;
    int write_allocate;
	int prefetch_num;
	int prefetch_stride;
	int workset_length;
	bool no_optim;
	int level;

	int ReplaceDecision(int set_idx);
	void ReplaceAlgorithm(uint64_t addr, int bytes, int read,
						char *content, int &hit, int &cycle, int prefetch_flag);
	void PrefetchAlgorithm(uint64_t addr, int hit_idx);

	Set *set;
	Storage *lower_;
	DISALLOW_COPY_AND_ASSIGN(Cache);
};

#endif //CACHE_CACHE_H_
