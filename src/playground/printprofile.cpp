#include <LocalParameters.h>
#include <Debug.h>
#include <DBReader.h>
#include <cstdlib>
#include <ctime>
#include <BaseMatrix.h>
#include <SubstitutionMatrix.h>

// include omp.h if def OPENMP
#ifdef OPENMP
#include <omp.h>
#endif

int printprofile(int argc, const char **argv, const Command &command)
{
    std::srand(std::time(nullptr));
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.rescoreMode = Parameters::RESCORE_MODE_ALIGNMENT;
    par.spacedKmer = false;
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    BaseMatrix *subMat = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, 0.0);

    std::vector<std::string> allProfileSeqs;
    allProfileSeqs.reserve(reader.getSize());

#pragma omp parallel num_threads(par.threads)
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        Sequence seq(par.maxSeqLen, reader.getDbtype(), subMat, 0, par.spacedKmer, false, false);

        std::vector<std::string> results;
#pragma omp for schedule(static)
        for (size_t i = 0; i < reader.getSize(); i++)
        {
            const char *seqData = reader.getData(i, thread_idx);
            seq.mapSequence(i, i, seqData, reader.getSeqLen(i));
            std::string realSeq;
            seq.extractProfileSequence(seqData, *subMat, realSeq);
            results.push_back(realSeq);
        }

#pragma omp critical
        {
            allProfileSeqs.insert(allProfileSeqs.end(), results.begin(), results.end());
        }
    }

    Debug(Debug::INFO) << "Randomly print profile sequences: \n";
    // randomly sample 20 sequences to print
    for (size_t i = 0; i < 20; i++)
    {
        size_t idx = std::rand() % allProfileSeqs.size();
        Debug(Debug::INFO) << allProfileSeqs[idx] << "\n";
    }

    return EXIT_SUCCESS;
}
