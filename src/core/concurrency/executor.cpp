#include "core/concurrency/executor.hpp"

#include <algorithm>
#include <mutex>
#include <thread>
#include <vector>

#ifdef PIRU_USE_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#endif

namespace piru::concurrency {
namespace {

class SerialTaskGroup : public TaskGroup {
public:
  void run(const std::function<void()>& fn) override { fn(); }
  void wait() override {}
};

class SerialExecutor : public Executor {
public:
  void parallel_for(std::size_t begin, std::size_t end, std::size_t grain,
                    const std::function<void(std::size_t)>& fn) override {
    const std::size_t step = grain == 0 ? 1 : grain;
    for (std::size_t i = begin; i < end; i += step) {
      const std::size_t limit = std::min(i + step, end);
      for (std::size_t j = i; j < limit; ++j) {
        fn(j);
      }
    }
  }

  int max_concurrency() const override { return 1; }
  std::string backend_name() const override { return "serial"; }
};

#ifdef PIRU_USE_TBB
class TbbTaskGroup : public TaskGroup {
public:
  void run(const std::function<void()>& fn) override { group_.run(fn); }
  void wait() override { group_.wait(); }

private:
  tbb::task_group group_;
};

class TbbExecutor : public Executor {
public:
  explicit TbbExecutor(int threads) : arena_(threads > 0 ? threads : tbb::task_arena::automatic) {}

  void parallel_for(std::size_t begin, std::size_t end, std::size_t grain,
                    const std::function<void(std::size_t)>& fn) override {
    const std::size_t step = grain == 0 ? 1 : grain;
    arena_.execute([&] {
      tbb::parallel_for(tbb::blocked_range<std::size_t>(begin, end, step),
                        [&](const tbb::blocked_range<std::size_t>& r) {
                          for (std::size_t i = r.begin(); i != r.end(); ++i) {
                            fn(i);
                          }
                        });
    });
  }

  int max_concurrency() const override { return arena_.max_concurrency(); }
  std::string backend_name() const override { return "oneTBB"; }

private:
  tbb::task_arena arena_;
};
#endif

}  // namespace

std::unique_ptr<Executor> make_executor(int threads) {
#ifdef PIRU_USE_TBB
  return std::make_unique<TbbExecutor>(threads);
#else
  return std::make_unique<SerialExecutor>();
#endif
}

std::unique_ptr<TaskGroup> make_task_group() {
#ifdef PIRU_USE_TBB
  return std::make_unique<TbbTaskGroup>();
#else
  return std::make_unique<SerialTaskGroup>();
#endif
}

}  // namespace piru::concurrency
