#include <BaseMatrix.h>
#include <Debug.h>
#include <Indexer.h>
#include <LocalParameters.h>
#include <SubstitutionMatrix.h>

int kmer2long(int argc, const char **argv, const Command &command) {
  LocalParameters &par = LocalParameters::getLocalInstance();
  par.parseParameters(argc, argv, command, true, 0,
                      LocalParameters::PARSE_VARIADIC);

  BaseMatrix *subMat = new SubstitutionMatrix(
      par.seedScoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2f);
  Indexer idx(subMat->alphabetSize - 1, 9);

  if (argc < 0) {
    Debug(Debug::ERROR)
        << "Too few arguments: please provide some inputs for conversion\n";
    return EXIT_FAILURE;
  }

  for (int i = 0; i < argc; i++) {
    Debug(Debug::INFO) << "Processing " << argv[i] << "\n";
    const unsigned char *kmer = (const unsigned char *)argv[i];
    Debug(Debug::INFO) << idx.int2index(kmer, 0, std::strlen(argv[i])) - 1
                       << "\n";
  }
  return EXIT_SUCCESS;
}
