#ifndef gem_BITMAP
#define gem_BITMAP

#define BITSHRINK(i) (i >> 6)
#define BITOFFSET(i) (i & 0x3f)

struct BitMap{
  int _size;
  unsigned long * _data;

  BitMap(int size, int type=0) { 
    int len = BITSHRINK(size)+1;
    this->_size = size;
    this->_data = new unsigned long [len];

    long INIT = (type==0 ? 0 : 0xffffffffffffffff);

    #pragma omp parallel for
    for (int i=0; i<len; i++) {
      this->_data[i] = INIT;
    }
    #pragma omp barrier
  }

  ~BitMap() {
    if(this->_data != NULL) {
      free(this->_data);
      this->_data=NULL;
    }
  }

  void clear(){
    size_t bm_size = BITSHRINK(this->_size);
    #pragma omp parallel for
    for (size_t i=0;i<=bm_size;i++) {
      this->_data[i] = 0;
    }
    #pragma omp barrier
  }

  void clear(int start, int count) {
    int end = start + count;
    int s_range = BITSHRINK(start);
    int s_offset= BITOFFSET(start);
    int e_range = BITSHRINK(end);
    int e_offset= BITOFFSET(end);
    if(s_range != e_range) {
      #pragma omp parallel for
      for (int i = s_range+1 ; i < e_range; i++) {
        this->_data[i] = 0;
      }
      #pragma omp barrier

      __sync_fetch_and_and(this->_data+s_range, ((1ul<<s_offset) - 1));
      __sync_fetch_and_and(this->_data+e_range, ~((1ul<<e_offset) - 1));
    } else {
      #pragma omp parallel for
      for (int i = start; i < end; i++) {
        clean_bit(i);
      }
      #pragma omp barrier
    }
  }

  bool get_bits(int start, int count) {
    int end = start + count;
    int s_range = BITSHRINK(start);
    int s_offset= BITOFFSET(start);
    int e_range = BITSHRINK(end);
    int e_offset= BITOFFSET(end);
    if (s_range != e_range) {
      for (int i = s_range + 1; i < e_range; i++) {
        if (this->_data[i] != 0) { return true; }
      }

      long s = this->_data[s_range];
      long e = this->_data[e_range];
      if ((s & ~((1ul<<s_offset) - 1)) > 0) return true;
      if ((e &  ((1ul<<e_offset) - 1)) > 0) return true; 
    } else {
      for (int i = start; i < end; i ++) {
        if (get_bit(i)) { return true; }
      }
    }
    return false;
  }

  bool get_bit(int i) {
    return bool(this->_data[BITSHRINK(i)] & (1ul << BITOFFSET(i)));
  }

  void set_bit(int i) {
    __sync_fetch_and_or(this->_data+BITSHRINK(i), 1ul<<BITOFFSET(i));
  }

  void clean_bit(int i) {
    __sync_fetch_and_and(this->_data+BITSHRINK(i), ~(1ul<<BITOFFSET(i)));
  }
};

#endif // gem_BITMAP
