#include <getopt.h>
#include <cassert>

#include "data_reader.h"
#include "realignment.h"
#include "group.h"

static int usage() {
	fprintf(stderr, "Usage: phase -b <BAM> -r <FASTA> -v <VCF> -o <Output>\n");
	fprintf(stderr, "  -b aligned reads in BAM format (indexed required)\n");
	fprintf(stderr, "  -r reference sequence for allele realignment in FASTA format (indexed required)\n");
	fprintf(stderr, "  -v heterozygous variants to phase in VCF format\n");
	fprintf(stderr, "  -o output file that phased results are written to (stdout)\n");
	fprintf(stderr, "  -c specify a chromosome to phase\n");
	return 1;
}

int main(int argc, char* argv[]) {
    if (argc == 1) return usage();    

    const char *bam_fn = nullptr, *vcf_fn = nullptr,  *ref_fn = nullptr, *output_fn = nullptr;
	const char *request_chromosome = nullptr;
	const char *resolution_fn = nullptr;
    int c;
    while ((c = getopt(argc, argv, "b:v:o:c:r:l:R:")) >= 0) {
		if (c == 'b') {
			bam_fn = optarg;
		} else if (c == 'v') {
			vcf_fn = optarg;
		} else if (c == 'o') {
			output_fn = optarg;
		} else if (c == 'c') {
			request_chromosome = optarg;
		} else if (c == 'r') {
			ref_fn = optarg;
		} else return usage();
	}

    if (bam_fn == nullptr) { fprintf(stderr, "ERR: please assign a BAM file of aligned reads\n"); return 1; }
	if (vcf_fn == nullptr) { fprintf(stderr, "ERR: please choose a VCF file to phase\n"); return 1; }
	if (ref_fn == nullptr) {
		fprintf(stderr, "ERR: please provide a reference sequence for detecting alleles by realignment\n" );
		return 1;
	}

    auto variant_table = input_vcf(vcf_fn, request_chromosome);

    FASTA_Reader ref_reader(ref_fn);
    
    // enumerate chrs
    for (int i = 0; i < variant_table.size; i ++) { // 遍历所有染色体	
        const auto &chr_name = variant_table.chromosomes[i];
		auto &snp_column = variant_table.variants[i]; // 对应染色体的所有snp
		fprintf(stderr, "Phase %ld SNPs on chromosome %s\n", snp_column.size(), chr_name.c_str());

		// Detecting alleles (include Realignment)
        int length = ref_reader.get_length(chr_name); assert(length > 0); // 获取染色体的长度
        char *sequence = ref_reader.get_contig(chr_name); // 获取染色体的序列
        std::vector<Read_Allele> read_row;
        read_row = detect_allele(bam_fn, chr_name, snp_column, length, sequence); // 检测 allele, 并进行realignment，返回的是所有 read 的 SNP 和 allele
        delete [] sequence;


		/**
		 * 思考一下要不要边分组边建树
		 * 这样的话，先有分组，再有建树
		 * 可以将树的根节点信息存储在组中
		 */
		// group SNPs by chunkL & chunkV
		auto groups = group_snps(snp_column);

		// create haplotype tree for each group
    }
}