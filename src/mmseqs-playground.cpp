#include "Command.h"
#include "DownloadDatabase.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char *binary_name = "mmseqs-playground";
const char *tool_name = "mmseqs-playground";
const char *tool_introduction = "MMseqs playground";
const char *main_author = "Matchy Lee";
const char *show_extended_help = "1";
const char *show_bash_info = "1";
const char *index_version_compatible = "16";

bool hide_base_commands = true;
bool hide_base_downloads = false;

void updateValidation();

void (*validatorUpdate)(void) = updateValidation;

std::vector<DatabaseDownload> externalDownloads = {};

LocalParameters &localPar = LocalParameters::getLocalInstance();

std::vector<struct Command> commands = {
        {"kmer2long",
                kmer2long,
                &localPar.kmer2long,
                COMMAND_MAIN,
                "",
                NULL,
                "Matchy Lee <matchy@snu.ac.kr>",
                "<i:kmer string>",
                CITATION_MMSEQS2,
                {}},
        {"long2kmer",
                long2kmer,
                &localPar.long2kmer,
                COMMAND_MAIN,
                "",
                NULL,
                "Matchy Lee <matchy@snu.ac.kr>",
                "<i:long integer>",
                CITATION_MMSEQS2,
                {}},
        {"printprofile",
                printprofile,
                &localPar.printprofile,
                COMMAND_MAIN,
                "",
                NULL,
                "Matchy Lee <matchy@snu.ac.kr>",
                "<i:profiledb> <i:seqDb>",
                CITATION_MMSEQS2,
                {{"profileDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::profileDb},
                        {"seqDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}}},
};

void updateValidation() {}
