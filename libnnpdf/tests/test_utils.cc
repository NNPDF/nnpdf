#include "catch.hpp"
#include "utils.h"

using namespace NNPDF;

TEST_CASE("Test joinpath", "[utils]"){
    REQUIRE(joinpath({""}) == "");
    REQUIRE(joinpath({}) == "");
    REQUIRE(joinpath({"a","b","cd"}) == "a/b/cd");
    REQUIRE(joinpath({"a","","b"}) == "a/b");
    REQUIRE(joinpath({"a","/b","c/"}) == "/b/c/");
    REQUIRE(joinpath({"/a","b","c/"}) == "/a/b/c/");
	auto comps = std::vector<std::string>{"a", "b", "c"};
	REQUIRE(joinpath(comps) == "a/b/c");

}
