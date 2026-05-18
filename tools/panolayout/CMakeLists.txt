# panolayout: standalone PG-SGD 1D layout tool.
# Header-only deps: third_party/dirtyzipf. Links pthread for std::async.

add_executable(panolayout panolayout.cpp)
target_include_directories(panolayout PRIVATE
  ${CMAKE_SOURCE_DIR}/third_party/dirtyzipf
)
target_compile_features(panolayout PRIVATE cxx_std_17)
find_package(Threads REQUIRED)
target_link_libraries(panolayout PRIVATE Threads::Threads)
