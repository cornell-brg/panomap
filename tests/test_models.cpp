#include <doctest/doctest.h>
#include <filesystem>

#include "io/models/model_factory.hpp"

TEST_CASE("built-in r9.4 model loads and exposes entries") {
  auto model = piru::io::load_builtin_model("r9.4");
  REQUIRE(model != nullptr);
  CHECK(model->k() == 6);

  double mean = 0.0;
  REQUIRE(model->lookup("AAAAAA", mean));
  CHECK(mean == doctest::Approx(86.486336));
}

TEST_CASE("built-in r10.4 model loads and exposes entries") {
  auto model = piru::io::load_builtin_model("r10.4");
  REQUIRE(model != nullptr);
  CHECK(model->k() == 9);

  double mean = 0.0;
  REQUIRE(model->lookup("AAAAAAAAA", mean));
  // Value from squigulator R10 model (raw pA)
  CHECK(mean == doctest::Approx(54.31598354));
}

TEST_CASE("load model from legacy r9-format file") {
  const auto src_dir = std::filesystem::path(__FILE__).parent_path();
  const auto path = src_dir / "data/models/r9_sample.model";
  auto model = piru::io::load_model_from_file(path.string());
  REQUIRE(model != nullptr);
  CHECK(model->k() == 6);
  double mean = 0.0;
  REQUIRE(model->lookup("AAAAGT", mean));
  CHECK(mean == doctest::Approx(82.866866));
}

TEST_CASE("load model from r10 two-column file") {
  const auto src_dir = std::filesystem::path(__FILE__).parent_path();
  const auto path = src_dir / "data/models/r10_sample.txt";
  auto model = piru::io::load_model_from_file(path.string());
  REQUIRE(model != nullptr);
  CHECK(model->k() == 9);
  double mean = 0.0;
  REQUIRE(model->lookup("AAAAAAACA", mean));
  CHECK(mean == doctest::Approx(-1.4318406));
}
