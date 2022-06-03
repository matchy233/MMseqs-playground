#include <BaseMatrix.h>
#include <DBReader.h>
#include <Debug.h>
#include <LocalParameters.h>
#include <SubstitutionMatrix.h>
#include <cstdlib>
#include <ctime>

// include omp.h if def OPENMP
#ifdef OPENMP

#include <omp.h>

#endif

int printprofile(int argc, const char **argv, const Command &command) {
    std::srand(std::time(nullptr));
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(
            par.db1.c_str(), par.db1Index.c_str(), par.threads,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> seqReader(
            par.db2.c_str(), par.db2Index.c_str(), par.threads,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    seqReader.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat = new SubstitutionMatrix(
            par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);

    std::vector<std::string> allProfileSeqs;
    allProfileSeqs.reserve(reader.getSize());

#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(par.maxSeqLen, reader.getDbtype(), subMat, 0, par.spacedKmer,
                     false, false);

        std::vector<std::string> results;
#pragma omp for schedule(static)
        for (size_t i = 0; i < reader.getSize(); i++) {
            const char *seqData = reader.getData(i, thread_idx);
            size_t queryKey = reader.getDbKey(i);
            size_t querySeqLen = reader.getSeqLen(i);
            seq.mapSequence(i, i, seqData, querySeqLen);
            std::string seqFromNumSeq;
            std::string seqFromWhileLoop;

            for (size_t i = 0; i < (querySeqLen) * Sequence::PROFILE_READIN_SIZE;
                 i += Sequence::PROFILE_READIN_SIZE) {
                seqFromWhileLoop.append(
                        1, subMat->num2aa[(int) seqData[i + Sequence::PROFILE_AA_SIZE]]);
            }

            for (size_t i = 0; i < querySeqLen; i++) {
                seqFromNumSeq.push_back(subMat->num2aa[seq.numSequence[i]]);
            }

            if (seqFromWhileLoop != seqFromNumSeq) {
                Debug(Debug::ERROR) << "ERROR: seqFromWhileLoop != seqFromNumSeq"
                                    << "\n";
                Debug(Debug::INFO) << "seqFromWhileLoop: " << seqFromWhileLoop << "\n";
                Debug(Debug::INFO) << "seqFromNumSeq: " << seqFromNumSeq << "\n";
            }

            // seq.extractProfileSequence(seqData, *subMat, realSeq);
            if (seqFromNumSeq.length() != querySeqLen) {
                Debug(Debug::ERROR)
                        << "Error: sequence length does not match for sequence " << i
                        << " with queryKey " << queryKey << ".\n";
                Debug(Debug::ERROR)
                        << "Expected equence length: " << querySeqLen << "\n";
                Debug(Debug::ERROR)
                        << "Actual sequence length: " << seqFromNumSeq.length() << "\n";
                EXIT(EXIT_FAILURE);
            }
            std::string seqFromSeqDb(seqReader.getData(i, thread_idx));
            if (seqFromNumSeq != seqFromSeqDb) {
                Debug(Debug::ERROR) << "ERROR: seqFromNumSeq != seqFromSeqDb"
                                    << "\n";
                Debug(Debug::INFO) << "seqFromSeqDb: " << seqFromSeqDb << "\n";
                Debug(Debug::INFO) << "seqFromNumSeq: " << seqFromNumSeq << "\n";
                EXIT(EXIT_FAILURE);
            }

            results.push_back(seqFromNumSeq);
        }

#pragma omp critical
        {
            allProfileSeqs.insert(allProfileSeqs.end(), results.begin(),
                                  results.end());
        }
    }

    Debug(Debug::INFO) << "Randomly print profile sequences: \n";
    // randomly sample 20 sequences to print
    for (size_t i = 0; i < 20; i++) {
        size_t idx = std::rand() % allProfileSeqs.size();
        Debug(Debug::INFO) << allProfileSeqs[idx] << "\n";
    }

    return EXIT_SUCCESS;
}
