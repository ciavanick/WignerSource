##
# @defgroup BashScripts Environment Setup Scripts
# @brief Bash scripts for setting up the Wigner simulation environment.
# @{
##
##
# @file wignerenv.sh
# @brief Sets environment variables for Wigner simulation and plotting.
#
# This script must be sourced manually in each new terminal session
# prior to running ROOT macros or the `wigneroot` executable.
#
# It updates:
# - `PATH` to include the `bin/` directory
# - `LD_LIBRARY_PATH` to include the `lib/` directory
#
# This script is already sourced automatically by `simulation.sh`,
# but must be sourced explicitly when running commands interactively.
#
# Example usage:
# @code
# source wignerenv.sh
# @endcode
##

export cpath=$(pwd)
export PATH="$PATH:$cpath/bin"
export LD_LIBRARY_PATH="$cpath/lib:${LD_LIBRARY_PATH:-}"
##
# @}
##