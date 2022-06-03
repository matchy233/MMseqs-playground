#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int kmer2long(int argc, const char **argv, const Command &command);

extern int long2kmer(int argc, const char **argv, const Command &command);

extern int printprofile(int argc, const char **argv, const Command &command);

extern int genkmer(int argc, const char **argv, const Command &command);

#endif
