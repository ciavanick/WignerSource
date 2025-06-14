#include "TSystem.h"
#include "TRint.h"

int main(int argc, char** argv) {
    // Set ROOT's internal search path to include our /wigner/lib
    gSystem->AddDynamicPath("./lib");  // for current dir
    gSystem->AddDynamicPath("../lib"); // for ../lib if executed from bin/
    gSystem->AddDynamicPath("/wigner/lib"); // absolute path in Docker

    // Preload dictionary
    gSystem->Load("libWignerUtils");

    TRint theApp("wigneroot", &argc, argv);
    theApp.Run();
    return 0;
}