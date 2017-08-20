#ifndef _BARRIER_
#define _BARRIER_
struct  x_barrier {
  volatile unsigned int count[2];
  volatile unsigned int sense;
  unsigned long expected;
// public:
};

  void init_barrier(struct x_barrier* sync, unsigned long expected_in)
  {
   sync->sense = 0;
  sync->expected = expected_in;
    sync->count[0] = 0;
    sync->count[1] = 0;
  }
  void wait_b(struct x_barrier* sync)
  {
    unsigned long sense_used = sync->sense;
    unsigned long arrived =
      __sync_fetch_and_add(&sync->count[sense_used], 1);
    if(arrived == (sync->expected - 1)) {
      sync->sense = 1 - sense_used; // Reverse sense
      sync->count[sense_used] = 0;
    }
    while(sync->count[sense_used] != 0) { NOP10();}
    __sync_synchronize(); // Also clobber memory
  }

#endif
