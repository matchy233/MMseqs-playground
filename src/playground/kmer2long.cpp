#include <Debug.h>
#include <LocalParameters.h>

int kmer2long(int argc, const char **argv, const Command &command) {
  LocalParameters &par = LocalParameters::getLocalInstance();
  par.parseParameters(argc, argv, command, true, 0,
                      LocalParameters::PARSE_VARIADIC);
  Debug(Debug::INFO) << par.db1 << "\n";
  return EXIT_SUCCESS;
}
