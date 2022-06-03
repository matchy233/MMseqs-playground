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

private:
    LocalParameters() : Parameters() {
        printprofile.push_back(&PARAM_THREADS);
        printprofile.push_back(&PARAM_K);
        printprofile.push_back(&PARAM_V);
    }

    LocalParameters(LocalParameters const &);

    ~LocalParameters() {};

    void operator=(LocalParameters const &);
};

#endif
