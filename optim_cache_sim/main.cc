#include "cache.h"
#include "memory.h"
#include "def.h"
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <cstring>
#include <stdint.h>

using namespace std;

double hit_rate1;
int main(int argc, char *argv[]) {

    int round = 0;
    char infile[50] = {0};
    int workset_length = 50;
    bool no_optim = false;
    int prefetch_num = 4;
    int prefetch_stride = 16;

    int opt;
    int option_index = 0;
    char optstring[20] = "i:r:l::p::s::";
    static struct option long_options[] = {
        {"no_optim", no_argument, NULL, 'n'},
        {0, 0, 0, 0}
    };
    
    while ( (opt = getopt_long(argc, argv, optstring, long_options, &option_index)) != -1)
    {
        switch(opt) {
            case 'i':
                strcpy(infile, optarg);
                cout << "Input file: " << infile << endl;
                break;
            case 'r':
                round = stoi(optarg);
                cout << "Number of iteration rounds:" << round << endl;
                break;
            case 'l':
                workset_length = stoi(optarg);
                break;
            case 'p':
                prefetch_num = stoi(optarg);
                break;
            case 's':
                prefetch_stride = stoi(optarg);
                break;
            case 'n':
                no_optim = true;
                cout << "Disable optimization metrics" << endl;
                break;
        }
    }

    cout << "Workset length: " << workset_length << endl;
    cout << "Prefetch lines: " << prefetch_num << endl;
    cout << "Prefetch Stride: " << prefetch_stride << endl;
    if (!no_optim) cout << "Enable optimization metrics" << endl;

    Memory m;
    Cache l1(L1_CACHE), l2(L2_CACHE);
    l1.SetConfig(l1_size, l1_asso, l1_block, l1_write_through, prefetch_stride,
                l1_write_allocate, prefetch_num, workset_length, no_optim);
    l2.SetConfig(l2_size, l2_asso, l2_block, l2_write_through, prefetch_stride,
                l2_write_allocate, prefetch_num, workset_length, no_optim);

	StorageLatency ml(m_hit_latency, m_bus_latency);
    m.SetLatency(ml);

    StorageLatency cl01(l1_hit_latency, l1_bus_latency);
    l1.SetLatency(cl01);
    StorageLatency cl02(l2_hit_latency, l2_bus_latency);
    l2.SetLatency(cl02);

    l1.SetLower(&l2);
    l2.SetLower(&m);

    StorageStats s;
    m.SetStats(s);
    l1.SetStats(s);
    l2.SetStats(s);

    int hit = 0, cycle = 0, trace_cnt = 0;
    char content[64];
    uint64_t addr;
    char method;
    char addr_hex[20] = {0};
    char Type;

    ifstream fin;
	for (int i = 0; i < round; i++) {
		fin.open(infile);
        while (fin >> Type >> addr_hex) {
            trace_cnt++;
            sscanf(addr_hex, "0x%llx", &addr);
            if (Type == 'w')
                l1.HandleRequest(addr, 0, FALSE, content, hit, cycle, FALSE);
            else
                l1.HandleRequest(addr, 0, TRUE, content, hit, cycle, FALSE);
            hit_rate1 = (double)l1.stats_.miss_cnt / (l1.stats_.access_cnt - l1.stats_.bypass_cnt);
        }
        fin.close();
    }

	double l1_missrate, l2_missrate, m_missrate;
    cout << endl;
	cout << "L1:"  << endl;
    l1.printStats(l1_missrate, false);
    cout << endl;
	cout << "L2:"  << endl;
    l2.printStats(l2_missrate, false);
    cout << endl;
	cout << "Main Memory:" << endl;
    m.printStats(m_missrate, true);

	cout << endl;
    cout << "Access Count: " << trace_cnt << endl;
    cout << "Total Access Time (in cycles): " << cycle << endl;
	double amat = l1_hit_latency + (l2_bus_latency + l2_hit_latency) * l1_missrate 
            + m_hit_latency * l1_missrate * l2_missrate;
	cout << "AMAT:" << amat << endl;
    cout << endl;
    return 0;
}
