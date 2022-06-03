#include <omp.h>
#include "LocalParameters.h"
#include "FileUtil.h"
#include "DBReader.h"
#include "SubstitutionMatrix.h"
#include "ExtendedSubstitutionMatrix.h"
#include "Indexer.h"
#include "KmerGenerator.h"

// define const unsigned char* as Kmer type
typedef size_t Kmer;

std::vector<Kmer> readAndGenerateKmer(std::string db, std::string dbIndex,
                                      LocalParameters &par, double similarKmerFactor = 1.5);


int genkmer(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, Parameters::PARSE_VARIADIC);

    std::vector<Kmer> kmerFromSeqDB = readAndGenerateKmer(par.db2, par.db2Index, par);
    std::vector<Kmer> kmerFromProfileDB = readAndGenerateKmer(par.db1, par.db1Index, par, 5);

    if (kmerFromProfileDB.size() < kmerFromSeqDB.size()) {
        Debug(Debug::WARNING) << "The number of kmer generated from profile database"
                              << "is fewer than that of sequence database!\n";
    }

    return EXIT_SUCCESS;
}

std::vector<Kmer> readAndGenerateKmer(std::string db,
                                      std::string dbIndex, LocalParameters &par,
                                      double similarKmerFactor) {

    int seqType = FileUtil::parseDbType(db.c_str());

    bool isProfile = Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_HMM_PROFILE);

    BaseMatrix *subMat;

    if (isProfile) {
        // the values are set according to lib/mmseqs/src/util/profile2seq.cpp
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.values.aminoacid().c_str(), 2.0f, 0.0);
    } else if (Parameters::isEqualDbtype(seqType, Parameters::DBTYPE_AMINO_ACIDS)) {
        subMat = new SubstitutionMatrix(par.seedScoringMatrixFile.values.aminoacid().c_str(), 8.0, -0.2f);
    } else {
        Debug(Debug::ERROR) << "Invalid input type (Support: amino acid, profile)\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "kmer threshold for profile: " << par.kmerScore.values.profile()
                       << "kmer threshold for sequence: " << par.kmerScore.values.sequence()
                       << "\n";

    DBReader<unsigned int> reader(db.c_str(), dbIndex.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    std::vector<Kmer> resultTable;

    size_t kmerCount = 0;
#pragma omp parallel for reduction(+:kmerCount) default(none) shared(reader, par)
    for (size_t i = 0; i < reader.getSize(); ++i) {
        size_t currentSequenceLength = reader.getSeqLen(i);
        //number of ungapped k-mers per sequence = seq.length-k-mer.size+1
        kmerCount += currentSequenceLength >= (size_t) par.kmerSize ? currentSequenceLength - par.kmerSize + 1 : 0;
    }

    Debug(Debug::INFO) << "Number of sequences: " << reader.getSize() << "\n";
    auto tableCapacity = static_cast<size_t>(similarKmerFactor * static_cast<double>(kmerCount + 1));
    resultTable.reserve(tableCapacity);

    int xIndex = subMat->aa2num[(int) 'X'];

    ScoreMatrix twoMatrix, threeMatrix;
    if (!isProfile) {
        twoMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 2);
        threeMatrix = ExtendedSubstitutionMatrix::calcScoreMatrix(*subMat, 3);
    }

    Debug::Progress progress(reader.getSize());
#pragma omp parallel default(none) \
shared(par, reader, subMat, progress, seqType, twoMatrix, threeMatrix, tableCapacity, resultTable, isProfile, xIndex)
    {
        unsigned int thread_idx = 0;
        unsigned int total_threads = 1;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
        total_threads = (unsigned int) omp_get_num_threads();
#endif
        Indexer idx(subMat->alphabetSize - 1, par.kmerSize);
        Sequence sequence(par.maxSeqLen, seqType, subMat, par.kmerSize, par.spacedKmer, false,
                          !isProfile);

        KmerGenerator kmerGenerator(par.kmerSize, subMat->alphabetSize - 1,
                                    isProfile ? par.kmerScore.values.profile() : par.kmerScore.values.sequence());

        if (isProfile && sequence.profile_matrix != nullptr) {
            kmerGenerator.setDivideStrategy(sequence.profile_matrix);
        } else {
            if (isProfile && sequence.profile_matrix == nullptr) {
                Debug(Debug::WARNING) << "The profile matrix retrieved from the profile database is empty!\n";
            }
            kmerGenerator.setDivideStrategy(&threeMatrix, &twoMatrix);
        }

        std::vector<Kmer> localTable;
        localTable.reserve(tableCapacity / total_threads);

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
//            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, (int) thread_idx);
            unsigned int seqLen = reader.getSeqLen(i);
            sequence.mapSequence(i, 0, data, seqLen);

            while (sequence.hasNextKmer()) {
                const unsigned char *kmer = sequence.nextKmer();

                int xCount = 0;
                for (int pos = 0; pos < par.kmerSize; ++pos) {
                    xCount += (kmer[pos] == xIndex);
                }

                if (xCount) {
                    continue;
                }

                if (isProfile) {
                    std::pair<size_t *, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
                    size_t lim = similarKmerList.second;
                    for (size_t j = 0; j < lim; ++j) {
                        localTable.emplace_back(similarKmerList.first[j]);
                    }
                    localTable.emplace_back(idx.int2index(kmer, 0, par.kmerSize));
                } else if (par.exactKmerMatching) {
                    localTable.emplace_back(idx.int2index(kmer, 0, par.kmerSize));
                } else {
                    // FIXME: too memory consuming when k = 11, need to adjust
                    //  (at least make the program does not terminate with an bad_alloc() error)
                    std::pair<size_t *, size_t> similarKmerList = kmerGenerator.generateKmerList(kmer);
                    size_t lim = similarKmerList.second;
                    for (size_t j = 0; j < lim; ++j) {
                        localTable.emplace_back(similarKmerList.first[j]);
                    }
                }
            }
        }

#pragma omp critical
        resultTable.insert(resultTable.end(), localTable.begin(), localTable.end());
    }
    return resultTable;
}
