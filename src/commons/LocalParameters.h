#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static void initInstance() { new LocalParameters; }

    static LocalParameters &getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters &>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter *> kmer2long;
    std::vector<MMseqsParameter *> long2kmer;
    std::vector<MMseqsParameter *> printprofile;
    std::vector<MMseqsParameter *> genkmer;

private:
    LocalParameters() : Parameters() {
        printprofile.push_back(&PARAM_THREADS);
        printprofile.push_back(&PARAM_K);
        printprofile.push_back(&PARAM_V);

        genkmer.push_back(&PARAM_EXACT_KMER_MATCHING);
        genkmer.push_back(&PARAM_SEED_SUB_MAT);
        genkmer.push_back(&PARAM_K);
        genkmer.push_back(&PARAM_K_SCORE);
        genkmer.push_back(&PARAM_SPACED_KMER_MODE);
        genkmer.push_back(&PARAM_MAX_SEQ_LEN);
        genkmer.push_back(&PARAM_COMPRESSED);
        genkmer.push_back(&PARAM_THREADS);
        genkmer.push_back(&PARAM_V);

        kmerSize = 9;
        kmerScore = 225;
    }

    LocalParameters(LocalParameters const &);

    ~LocalParameters() {};

    void operator=(LocalParameters const &);
};

#endif
