#include <sstream>
#include <algorithm>
#include <memory.h>
#include <cassert>
#include <fstream>

#include "data_reader.h"
#include "htslib/sam.h"

void VCF_Header::add_line(const char *buf) {
    if(buf[0] == '#' and buf[1] == '#'){
        std::string key;
        std::string value;
        int i = 2;
        for(; buf[i] != '='; i++) {
            key += buf[i];
        }
        i++;
        for(; buf[i] != '\n'; i++) {
            value += buf[i];
        }
        line_n ++;
        keys.push_back(key);
        values.push_back(value);
    } else if(buf[0] == '#') {
        auto s = split_str(buf, '\t');
        if (s.size() - VCF_SAMPLE != 1) {
			fprintf(stderr, "ERR: %ld samples are found, but KSNP is designed for a single individual\n",
		        s.size() - VCF_SAMPLE);
		}
        sample_line = buf;
    }
}

void VCF_Header::annotate_for_phasing() {
	bool has_gt = false, has_ps = false;
	for (int i = 0; i < line_n; i++) {
		if (keys[i] == "FORMAT") {
			if (values[i].find("ID=GT") != std::string::npos) has_gt = true;
			if (values[i].find("ID=PS") != std::string::npos) has_ps = true;
		}
	}
	if (not has_gt) add_line("<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	if (not has_ps) add_line("<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">");
}

std::vector<std::string> VCF_Header::to_str() {
	std::vector<std::string> ret;
	for (int i = 0; i < line_n; i++) {
		ret.push_back("##" + keys[i] + "=" + values[i]);
	}
	ret.push_back(sample_line);
	return ret;
}

Variant_Table input_vcf(const char *fn, const char* chromosome) {
    gzFile in = gzopen(fn, "r");
    if(in == nullptr) {
        fprintf(stderr, "ERR: open vcf file %s falied.\n", fn);
        std::abort();
    }

    Variant_Table vt; vt.size = 0;
    const int SNP_BUF_SIZE = 4 * 1024 * 1024;
    char *buf = new char[SNP_BUF_SIZE];
    std::map<std::string, int> dict;
    
    while(gzgets(in, buf, SNP_BUF_SIZE) != nullptr) {
    //     // 首先读header, VCF_Header的add_line处理的字符串末尾是\0，
    //     // 但是gzgets会读到\n并保存，所以首先将读到的一行末尾改成\n.
        if(buf[0] == '#') { vt.header.add_line(buf); continue;}
        
        // 录入染色体
        auto fields = split_str(buf, '\t');
        const auto &chr = fields[VCF_CHROM];
        if (chromosome and chr != std::string(chromosome)) continue;

        int pos = stoi(fields[VCF_POS]);

        if(fields[VCF_REF].size() != 1 or fields[VCF_ALT].size() != 1) { continue; } // not a SNV
        char ref = fields[VCF_REF][0];
        char alt = fields[VCF_ALT][0];

        // 因为在Variant_Table中,chr和SNPs是以index作对应的,所以在这里进行中间情况的保存
        if(dict.find(chr) == dict.end()) { // 当前chr第一次出现
            vt.chromosomes.push_back(chr);
            vt.variants.emplace_back(std::vector<SNP>());
            vt.size++;
            dict[chr] = vt.size - 1;
        }
        
        vt.variants[dict[chr]].emplace_back(SNP(pos, ref, alt, buf));

    }

    delete[] buf;
    gzclose(in);

    if (vt.size == 0) {
		fprintf(stderr, "ERR: input no variants to phase\n");
		std::abort();
	}

	int snp_n = 0;
	for (auto &variant: vt.variants) {
		snp_n += variant.size();
		std::sort(variant.begin(), variant.end());
	}
    // std::cout << snp_n << std::endl;

    return vt;
}

std::vector<std::string> split_str(const char *s, char sep) {
    std::vector<std::string> ret;
    std::string temp;
    for(int i = 0; true; i++) {
        if(s[i] == '\0') {
            if(!temp.empty()) ret.push_back(temp);
            break;
        } else if(s[i] == sep) {
            if(!temp.empty()) ret.push_back(temp);
            temp = "";
            continue;
        }
        temp += s[i];
    }

    return ret;
}

FASTA_Reader::FASTA_Reader(const std::string &fn) {
    ref_fp = gzopen(fn.c_str(), "r");
    if(ref_fp == nullptr) {
        fprintf(stderr, "ERR: open fasta file %s failed. \n", fn.c_str());
        std::abort();
    }

    std::ifstream in(fn + ".fai");
    if(!in.is_open()) {
        fprintf(stderr, "ERR: reference sequence is not indexed. Use `samtools faidx`\n");
		std::abort();
    }

    auto *line_buf = new char[4 * 1024];
    while ( in.getline(line_buf, 4 * 1024)) {
        auto idxs = split_str(line_buf, '\t');
        std::string chr = idxs[0];

        FAIDX_Contig t;
        int length = std::stoi(idxs[1]);
        int offset = std::stoi(idxs[2]);
        int line_bases = std::stoi(idxs[3]);
        int line_width = std::stoi(idxs[4]);
        
        t.length = length; t.offset = offset; t.line_bases = line_bases; t.line_width = line_width;

        dict[chr] = t;
    }

    delete [] line_buf;
    in.close();
}


void newlinein(const char *str) {
    for(int i = 0; ;i++) {
        if(str[i] == '\n') {
            fprintf(stderr, "new line!\n");
            std::abort();
        }else if(str[i] == '\0') {
            fprintf(stderr, "\\0 \n");
            std::abort();
        }else {
            std::cout << str;
        }
    }
}

char * FASTA_Reader::get_contig(const std::string &chr_name) {

    const auto &faidx = dict[chr_name];
    char * ref = new char[faidx.length + 5]; int length = 0;
    gzseek(ref_fp, faidx.offset, SEEK_SET);
    char * buffer = new char[faidx.line_width + 5];

    // 不能用类似gzgets.. != nullptr.
    // 因为我们只要目标染色体的序列
    // 并且我们知道其偏移量和碱基长度,所以可以用长度来控制序列的读取
    while (length < faidx.length) {
        gzgets(ref_fp, buffer, faidx.line_width + 5);
        int n = std::min(faidx.line_bases, faidx.length - length);

        memcpy(ref + length, buffer, n);
        length += n;
    }

    delete[] buffer;
    assert(length == faidx.length);
    return ref;
}
