#include "data_reader.h"
#include <fstream>

int main(int argc, char **argv) {
    const char * fn = "./data/variants.vcf";
    const char * fa = "./data/ref.fa";

    // gzFile f = gzopen(fn, "r");

    auto table = input_vcf(fn, nullptr);
    auto fasta = FASTA_Reader(fa);




    std::ofstream of("out.txt", std::ios::out);
    char *sequence = fasta.get_contig("1");
    of << sequence;


    /**
     * TEST VARIANT_TABLE rel
     */
    // for(const auto &snps : table.variants) {
    //     for(const auto &snp : snps) {
    //          std::cout << snp.line << std::endl;
    //     }
    // }
    // auto ss = table.header.to_str();

    // for(const auto &s : ss) {
    //     std::cout << s << std::endl;
    // }

    
    // for(auto it = fasta.dict.begin(); it != fasta.dict.end(); ++it) {
    //     std::cout << "Key: " << it->first << ", Value: " << it->second.to_str() << std::endl;
    // }

    return 0;
}