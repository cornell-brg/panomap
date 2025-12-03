#include <doctest/doctest.h>

#include <string>

#include "io/reads/read_provider_factory.hpp"
#include "io/reads/read_provider.hpp"

TEST_CASE("slow5 provider reads bundled blow5") {
#ifdef PIRU_HAS_SLOW5
    const std::string path = "tests/data/HLA/test_reads/quick_r9_2k.blow5";
    auto provider = piru::io::make_read_provider(path);
    REQUIRE(provider != nullptr);
    CHECK(provider->get_format_name() == "slow5");

    piru::io::RawRead read;
    std::size_t count = 0;
    while (provider->get_next(read)) {
        ++count;
        CHECK_FALSE(read.read_id.empty());
        CHECK(read.len_raw_signal > 0);
    }
    CHECK(count == 5);
#else
    MESSAGE("PIRU_HAS_SLOW5 not set; skipping slow5 read test");
    CHECK(true);
#endif
}

// Simple mock provider to test interface usage without slow5.
class MockProvider : public piru::io::ReadProvider {
public:
    bool get_next(piru::io::RawRead& read) override {
        if (done_) return false;
        read.read_id = "mock";
        read.raw_signal = {1, 2, 3};
        read.len_raw_signal = read.raw_signal.size();
        done_ = true;
        return true;
    }
    void reset() override { done_ = false; }
    std::string get_format_name() const override { return "mock"; }

private:
    bool done_{false};
};

TEST_CASE("mock read provider") {
    MockProvider p;
    piru::io::RawRead read;
    CHECK(p.get_next(read));
    CHECK(read.len_raw_signal == 3);
    CHECK_FALSE(p.get_next(read));
    p.reset();
    CHECK(p.get_next(read));
}
