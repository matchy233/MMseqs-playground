#include <LocalParameters.h>
#include <Debug.h>

int printprofile(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0,
                        LocalParameters::PARSE_VARIADIC);


}
