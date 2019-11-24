#include "cache.h"
#include "memory.h"
#include "def.h"
#include <fstream>
#include <cstring>
#include <getopt.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

    Memory m;
    Cache l1;
    l1.SetLower(&m);

    char infile[20];
    int size = 0;
    int associativity = 0;
    int block_size = 0;
    int hit_latency = 0;
    int set_num = 0;
    int write_through = FALSE;
    int write_allocate = FALSE;
    int debug = FALSE; // dubug mode flag

    StorageStats s;
    s.access_time = 0;
    m.SetStats(s);
    l1.SetStats(s);
    
    //default latency
    StorageLatency ml;
    ml.hit_latency = 100;
    m.SetLatency(ml);

    StorageLatency ll;
    ll.hit_latency = 10;
    l1.SetLatency(ll);

    int opt;
    int option_index = 0;
    char optstring[20] = "i:c:s:b:l:";
    static struct option long_options[] = {
        {"write_through", no_argument, NULL, 'w'},
        {"write_allocate",  no_argument, NULL, 'a'},
        {"debug", no_argument, NULL, 'd'},
        {0, 0, 0, 0}
    };
    
    while ( (opt = getopt_long(argc, argv, optstring, long_options, &option_index)) != -1)
    {
        switch(opt) {
            case 'i':
                strcpy(infile, optarg);
                cout << "Input file: " << infile << endl;
                break;
            case 'c':
                size = stoi(optarg);
                cout << "Size of cache (in KB):" << size << endl;
                break;
            case 's':
                associativity = stoi(optarg);
                cout << "Associativity: " << associativity << endl;
                break;
            case 'b':
                block_size = stoi(optarg);
                cout << "Block size (in Byte): " << block_size << endl;
                break;
            case 'l':
                hit_latency = stoi(optarg);
                cout << "Cache hit latency (in Cycles): " << hit_latency << endl;
                break;
            case 'w':
                write_through = TRUE;
                cout << "Write through enabled" << endl;
                break;
            case 'a':
                write_allocate = TRUE;
                cout << "Write allocate enabled" << endl;
                break;
            case 'd':
                debug = TRUE;
                cout << "Debug mode enabled" << endl;
        }
    }

    if (!write_through) cout << "Write back enabled" << endl;
    if (!write_allocate) cout << "Write Allocate disabled" << endl;
    if (!debug) cout << "Debug mode disabled" << endl;

	ll.hit_latency = hit_latency;
    //initialize
    size = 1024 * size; // turn KBs to Bytes
    set_num = size / (associativity * block_size); // number of sets
    cout << "Number of Sets: " << set_num << endl;
    l1.SetConfig(size, associativity, set_num, write_through, write_allocate, debug);

    int hit = 0, cycle = 0, trace_cnt = 0, replace_cnt = 0;
    char content[64] = {0};
    uint64_t addr;
    char addr_hex[20] = {0};
    char Type;

    // read the traces
    ifstream fin;
    fin.open(infile);
    while (fin >> Type >> addr_hex) {
        trace_cnt++;
        sscanf(addr_hex, "0x%llx", &addr);
        if (debug)
            printf("Query: %c 0x%llx\n", Type, addr);
        if (Type == 'w')
            l1.HandleRequest(addr, 0, FALSE, content, hit, cycle, replace_cnt);
        else
            l1.HandleRequest(addr, 0, TRUE, content, hit, cycle, replace_cnt);
    }
    
    cout << endl;
    cout << "Total Access: " << trace_cnt << endl;
    cout << "Cache Hit: " << hit << endl;
    cout << "Hit Rate: " << ((double)hit / (double)trace_cnt) << endl;
    cout << "Miss Rate: " << 1 - ((double)hit / (double)trace_cnt) << endl;
    cout << "Total Replacement: " << replace_cnt << endl;
    l1.GetStats(s);
    cout << "Cache access time (in Cycles): " << s.access_time << endl; 
    cout << "Cache access time (in ms): " << (double)s.access_time / cpu_freq * 1000 << endl;
    m.GetStats(s);
    cout << "Memory access time (in Cycles): " << s.access_time << endl;
    cout << "Memory access time (in ms): " << (double)s.access_time / cpu_freq * 1000 << endl;
    cout << "Total access time (in Cycles): " << cycle << endl;
    cout << "Total access time (in ms): " << (double)cycle / cpu_freq * 1000 << endl;
    fin.close();
    return 0;
}
