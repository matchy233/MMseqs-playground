#include <BaseMatrix.h>
#include <Debug.h>
#include <Indexer.h>
#include <LocalParameters.h>
#include <SubstitutionMatrix.h>
#include <Util.h>

int long2kmer(int argc, const char **argv, const Command &command) {
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
    size_t kmeraAsLong = Util::fast_atoi<size_t>(argv[i]);
    Debug(Debug::INFO) << "kmeraAsLong: " << kmeraAsLong << "\n";
    idx.printKmer(kmeraAsLong, 9, subMat->num2aa);
    Debug(Debug::INFO) << "\n";
  }
  return EXIT_SUCCESS;
}
