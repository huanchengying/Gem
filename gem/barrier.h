#ifndef GEM_BARRIER
#define GEM_BARRIER

#include <thread>
#include <iostream>
#include <atomic>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <unistd.h>

namespace gem
{

class CVBarrier {

public:
  CVBarrier(const CVBarrier&) = delete;
  CVBarrier& operator=(const CVBarrier&) = delete;
  explicit CVBarrier(unsigned int count) 
  : _count(count), _generation(0), _count_reset_value(count)
  {}
  
  void Wait() {
    std::unique_lock< std::mutex > lock(_mutex);
    unsigned int gen = _generation;
    if (--_count == 0){
      _generation++;
      _count = _count_reset_value;
      _cond.notify_all();
      return;
    }
    while (gen == _generation) _cond.wait(lock);
  }

private:
  std::mutex _mutex;
  std::condition_variable _cond;
  unsigned int _count;
  unsigned int _generation;
  unsigned int _count_reset_value;
};

} // namespace gem
#endif // WONDRLAND_BARRIER