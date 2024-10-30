#ifndef GROUP_H
#define GROUP_H

#include <vector>
#include "data_reader.h"

class Group {
public:
    static const int MAX_CHUNK_LENGTH = 100000; // chunkL
    static const int MAX_CHUNK_VARIANTS = 100; // chunkV
    std::vector<SNP> snp_list; // 存储 SNP 的列表
    int start; // 第一个 SNP 的位置
    int end; // 最后一个 SNP 的位置
    int size; // SNP 的数量


};
typedef std::vector<Group> Groups;

// 接收染色体上所有 SNP，返回所有 SNP 的组
Groups group_snps(const std::vector<SNP> &snps);

#endif

