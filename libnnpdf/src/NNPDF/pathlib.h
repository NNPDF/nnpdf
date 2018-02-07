#pragma once
#include <string>

/**
 * pthlib.h
 *
 * Locate the nnprofile configuration file and extract basic
 * locations needed for other C++ projecs.
 */

namespace NNPDF
{
/**
 * @brief get_profile_path.
 * This is detemined by the environment variable
 * NNPDF_PROFILE_PATH if present, and otherwise by the compile time
 * constant DEFAULT_NNPDF_PROFILE_PATH.
 * @return The location of the NNPDF configuration file.
 */
std::string get_profile_path();

/**
 * @brief get_data_path
 * @return The location of the nnpdfcpp data path,
 * containing experimental data and theory settings.
 */
std::string get_data_path();

/**
 * @brief get_results_path
 * @return The location of the nnpdfcpp results path, containing
 * completed fits.
 */
std::string get_results_path();

/**
 * @brief get_config_path
 * @return The location of the nnpdfcpp config path, containing
 * default configurations for plotting and fiatlux
 */
std::string get_config_path();
} // namespace NNPDF
