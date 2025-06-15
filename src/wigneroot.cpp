/**
 * @defgroup WignerApp ROOT Interactive Entry Point
 * @brief Code for launching the ROOT interactive session with Wigner classes preloaded.
 *
 * This file defines the main function that starts a ROOT session using TRint,
 * sets up the library search paths, and dynamically loads the `libWignerUtils` shared library.
 * It is intended as the entry point for interactive macro usage (e.g., wignersim, makeplots).
 * @{
 */

#include "TSystem.h"
#include "TRint.h"

/**
 * @file wigneroot.cpp
 * @brief Entry point for launching the interactive Wigner analysis environment with ROOT.
 *
 * This executable initializes the ROOT interpreter, sets up custom library search paths,
 * and loads the WignerUtils shared library (libWignerUtils.so), which contains tools for
 * Wigner function evaluation and deuteron coalescence probability calculations.
 *
 * Before running this executable, make sure to set the correct environment using the
 * provided script `wignerenv.sh`, which configures `PATH` and `LD_LIBRARY_PATH`.
 *
 * Example usage:
 * @code
 *   source wignerenv.sh
 *   wigneroot
 * @endcode
 *
 * This file is built via the CMake target `wigneroot` and links to the shared library
 * `WignerUtils`, which is built from CWignerSource.cpp and CWignerUtils.cpp.
 */
 
/**
 * @brief Main function to launch ROOT interpreter with WignerUtils and WignerSource loaded.
 *
 * Configures dynamic library paths, loads the WignerUtils shared library,
 * and starts an interactive ROOT session via TRint.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit code.
 */
int main(int argc, char** argv) {
    // Add directories where ROOT can search for shared libraries
    gSystem->AddDynamicPath("./lib");        ///< Local lib path for development
    gSystem->AddDynamicPath("../lib");       ///< Relative lib path if run from bin/
    gSystem->AddDynamicPath("/wigner/lib");  ///< Absolute path (e.g., in Docker)

    // Load the WignerUtils shared library with ROOT dictionary
    gSystem->Load("libWignerUtils");

    // Launch interactive ROOT session
    TRint theApp("wigneroot", &argc, argv);
    theApp.Run();

    return 0;
}
/// @}