#ifndef READER_H
#define READER_H

#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <zlib.h>
#include <cstring>

const int VCF_CHROM  = 0;
const int VCF_POS    = 1;
const int VCF_ID     = 2;
const int VCF_REF    = 3;
const int VCF_ALT    = 4;
const int VCF_QUAL   = 5;
const int VCF_FILTER = 6;
const int VCF_INFO   = 7;
const int VCF_FORMAT = 8;
const int VCF_SAMPLE = 9;

struct VCF_Header {
private:
    int line_n;
    std::vector<std::string> keys; /** 例如INFO,FILTER,FORMAT这些 */
    std::vector<std::string> values; /** {keys}=后面的内容 */
    std::string sample_line;

public:
    VCF_Header(): line_n(0) {};

    /** Parse a header line in VCF file. shoule extract key & value */
    void add_line(const char* buf);
    void annotate_for_phasing();
    std::vector<std::string> to_str();

};

struct SNP {
    int pos;
    char ref; // 该位点在ref上对应base
    char alt; // vcf对应base
    char *line;

    // 记录覆盖该SNP的read以及read上该位点的base与ref/alt一样，或者没有。
    std::vector<int> rid; /** Read indices on this SNP */ // 因为read会被存储下来
    std::vector<int> allele; /** Alleles of reads. 0:REF 1:ALT -1:gap */

    /** Phased result */
    int ps;
    int gt; /** GT=0 for 0|1; GT=1 for 1|0; GT=-1 for unknown */

    SNP(int p, char r, char a, const char* l): pos(p), ref(r), alt(a) {
        ps = -1;
        gt = -1;
        if(l) line = strdup(l);
        else line = nullptr;
    }

    inline void add_read(int r, int a) {
        rid.push_back(r);
        allele.push_back(a);
    }

    inline int size() const { return (int)rid.size(); }
    
    bool operator < (const SNP &o) const { return this->pos < o.pos; }
};

struct Variant_Table {
    int size;
    VCF_Header header;
    std::vector<std::string> chromosomes;
    std::vector<std::vector<SNP>> variants; // 对应各chr上的SNPs。
};

Variant_Table input_vcf(const char *fn, const char *chromosome);

std::vector<std::string> split_str(const char *s, char sep);

/**************
 *   FASTA    *
 **************/
struct FAIDX_Contig {
    int length;
    int offset;
    int line_bases;
    int line_width;

    FAIDX_Contig(): length(0), offset(0), line_bases(0), line_width(0) {};

    inline std::string to_str() {
        std::ostringstream oss;
        oss << length << ',' << offset << ',' << line_bases << ',' << line_width;
        return oss.str();
    }
};

class FASTA_Reader {
public:
    gzFile ref_fp;
    std::map<std::string, FAIDX_Contig> dict; /** Dictionary from chromosome to offset */

public:
    explicit FASTA_Reader(const std::string &fn);

    inline int get_length(const std::string &chr_name) {
        return dict.find(chr_name) != dict.end() ? dict[chr_name].length : -1;
    }

    char *get_contig(const std::string &chr_name);

    inline void close() { gzclose(ref_fp); }
};


/**************
 *    BAM     *
 **************/
struct Allele_Call {
    /** 
     * 一个等位基因（附加存在），根据后面的用途：
     * SNPs会分别按reads和vcf文件进行存储，即记录各read上的SNPs信息 + 保存vcf文件中SNPs的信息
     * 因此该数据结构作为一个附加的存在，需记录
     * 1. 其在query read上的位置，也就是read上第几个snp
     * 2. 在vcf文件SNPs上的位置，同样，也是基于index保存
     * 3. 该等位基因信息，方便起见， 0表示与REF相同，1表示与read支持的ALT相同，-1表示gap/unknown
    */
   int que_pos; //position on query read
   int snp_idx; //position on vcf snps
   int allele;
   Allele_Call(int q, int s, int a): que_pos(q), snp_idx(s), allele(a) {}
};
typedef std::vector<Allele_Call> Read_Allele; // 以read为单位保存SNP的形式

/**
 * Detect alleles by realignment
 * 这个函数会进行realignment。
 * @param bam_fn aligned reads
 * @param chr_name chromosome name
 * @param snps variants to phase
 * @param len reference sequence length
 * @param seq reference sequence
 * @return all alleles on informative reads
 */
std::vector<Read_Allele> detect_allele(const char *bam_fn, const std::string &chr_name,
									   std::vector<SNP> &snps, int len, const char *seq);

#endif

