#ifndef BITMAP_HPP
#define BITMAP_HPP

#define WORD_OFFSET(i) ((i) >> 6)
#define BIT_OFFSET(i) ((i) & 0x3f)
template <class T>
inline bool cas(T * ptr, T old_val, T new_val) {
  if (sizeof(T) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&old_val), *((long*)&new_val));
  } else if (sizeof(T) == 4) {
    return __sync_bool_compare_and_swap((uint32_t*)ptr, *((uint32_t*)&old_val), *((uint32_t*)&new_val));
  } else {
    assert(false);
  }
}
class Bitmap {
public:
  size_t size;
  unsigned long * data;
  Bitmap() : size(0), data(NULL) { }
  Bitmap(size_t size) : size(size) {
    data = new unsigned long [WORD_OFFSET(size)+1];
    clear();
  }
  ~Bitmap() {
    delete [] data;
  }
  void clear() {
    size_t bm_size = WORD_OFFSET(size);
    parallel_for (size_t i=0;i<=bm_size;i++) {
      data[i] = 0;
    }
  }
  void fill() {
    size_t bm_size = WORD_OFFSET(size);
    parallel_for (size_t i=0;i<bm_size;i++) {
      data[i] = 0xffffffffffffffff;
    }
    data[bm_size] = 0;
    for (size_t i=(bm_size<<6);i<size;i++) {
      data[bm_size] |= 1ul << BIT_OFFSET(i);
    }
  }
  unsigned long get_bit(size_t i) {
    return data[WORD_OFFSET(i)] & (1ul<<BIT_OFFSET(i));
  }
  void set_bit(size_t i) {
    __sync_fetch_and_or(data+WORD_OFFSET(i), 1ul<<BIT_OFFSET(i));
  }
};

#endif

