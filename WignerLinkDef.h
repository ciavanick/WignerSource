/** 
 * \file WignerLinkDef.h
 * \brief ROOT dictionary configuration for wignerSource and wignerUtils.
 * 
 * This file contains `#pragma link` directives used by ROOT's CINT/Cling 
 * to generate class dictionaries for I/O and introspection.
 */

 #ifdef __CLING__
 #pragma link off all globals;   ///< Disable global variable linking
 #pragma link off all classes;   ///< Disable class linking by default
 #pragma link off all functions; ///< Disable function linking by default
 
 #pragma link C++ class wignerSource+; ///< Enable ROOT dictionary for wignerSource
 #pragma link C++ class wignerUtils+;  ///< Enable ROOT dictionary for wignerUtils
 #endif