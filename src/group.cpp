#include "group.h"

Groups group_snps(const std::vector<SNP> &snps) {
    Groups groups;
    const auto chunkL = Group::MAX_CHUNK_LENGTH;
    const auto chunkV = Group::MAX_CHUNK_VARIANTS;

    int i = 0;
    while (i < snps.size()) {
        Group group;
        for (int j = 0; j < chunkV && i + j < snps.size(); j++) {
            group.snp_list.push_back(snps[i + j]);
        }
        groups.push_back(group);
        i += chunkV;
    }
}

