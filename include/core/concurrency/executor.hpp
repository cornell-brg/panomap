#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <string>

namespace panomap::concurrency {

// Lightweight abstraction over the concurrency backend (oneTBB-backed by default).
class TaskGroup {
public:
  virtual ~TaskGroup() noexcept(false) = default;
  virtual void run(const std::function<void()>& fn) = 0;
  virtual void wait() = 0;
};

class Executor {
public:
  virtual ~Executor() noexcept(false) = default;
  virtual void parallel_for(std::size_t begin, std::size_t end, std::size_t grain,
                            const std::function<void(std::size_t)>& fn) = 0;
  virtual int max_concurrency() const = 0;
  virtual std::string backend_name() const = 0;
};

// Factory helpers (backend chosen at compile-time).
std::unique_ptr<Executor> make_executor(int threads = -1);
std::unique_ptr<TaskGroup> make_task_group();

}  // namespace panomap::concurrency
