#ifndef println

#define println(...) { if(me == 0) { printf(__VA_ARGS__); printf("\n"); } }
#define println_all(...) { printf(__VA_ARGS__); printf("\n"); }
#define print(...) { printf(__VA_ARGS__); }

#endif
