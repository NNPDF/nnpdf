#include "NNPDF/exceptions.h"
#include "NNPDF/pathlib.h"
#include <cstdlib>
#include <yaml-cpp/yaml.h>

namespace NNPDF
{

namespace
{

const YAML::Node &load()
{
    auto p = get_profile_path();
    try {
        const static auto res = YAML::LoadFile(p);
        return res;
    } catch (YAML::Exception &e) {
        throw FileError("pathlib", "Could not process profile file '" + p
                                       + "': " + e.what());
    }
}

template <typename T> T get_key(std::string param)
{
    auto node = load()[param];
    // yaml-cpp error messages when a key is not found are
    // pretty horrible, so we write our own.
    if (!node) {
        throw RuntimeException(
            "pathlib", "'" + param + "' " + "key not found in profile file '"
                           + get_profile_path() + "'");
    }
    return node.as<T>();
}

} // namespace

std::string get_profile_path()
{
    char *env_path = std::getenv("NNPDF_PROFILE_PATH");
    if (env_path) {
        return env_path;
    }
    // This is defined at compile time
    return DEFAULT_NNPDF_PROFILE_PATH;
}

std::string get_data_path() { return get_key<std::string>("data_path"); }

std::string get_results_path() { return get_key<std::string>("results_path"); }

std::string get_config_path() { return get_key<std::string>("config_path"); }

} // namespace NNPDF

