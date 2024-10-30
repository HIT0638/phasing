#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, getopt
import pandas as pd
import numpy as np
import subprocess

columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]


def main(argv):
    vcf_file = ''
    try:
        # 修改选项解析格式，允许 -v 后带有参数
        opts, args = getopt.getopt(argv, "v:", ["ifile="])
    except getopt.GetoptError:
        print('filter.py -v <vcffile>')
        sys.exit(2)

    # 解析选项和参数
    for opt, arg in opts:
        if opt == '-v':
            vcf_file = arg

    vcf_data = pd.read_csv(vcf_file, comment='#', sep='\t', header=None)
    vcf_data.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

    
    # filter_by_qual(vcf_data)
    # filter_by_dp(vcf_data)
    filter_by_1kgp3(vcf_data)



def filter_by_qual(vcf_data):
    # 按QUAL排序，计算minQUAL和maxQUAL
    qual_values = vcf_data["QUAL"].sort_values() # 获取所有QUAL值
    min_qual, max_qual = qual_values.min(), qual_values.max() # 获取最小和最大的QUAL值

    # 将[minQUAL, maxQUAL]均分为20个区间
    intervals = np.linspace(min_qual, max_qual, 21) # 将区间均分为20个      
    # 找到30%位数对应的区间位置
    threshold = np.percentile(qual_values, 30) # 找到30%位数对应的QUAL值
    interval_index = np.searchsorted(intervals, threshold) # 找到对应的区间位置

    # 过滤掉第1到i-1个区间
    qual_filtered_data = vcf_data[vcf_data["QUAL"] >= intervals[interval_index - 1]] # 过滤掉第1到i-1个区间

    # 保存过滤后的VCF文件
    qual_filtered_data.to_csv("filtered_by_qual.vcf", sep='\t', index=False)

def filter_by_dp(vcf_data):
    # 从INFO字段提取DP值（假设INFO字段格式类似"DP=100;..."）
    vcf_data["DP"] = vcf_data["INFO"].str.extract(r"DP=(\d+)").astype(float)  

    dp_values = vcf_data["DP"].sort_values()
    min_qual, max_qual = dp_values.min(), dp_values.max()

    intervals = np.linspace(min_qual, max_qual, 21)
    threshold = np.percentile(dp_values, 30)
    interval_indx = np.searchsorted(intervals, threshold)

    dp_filtered_data = vcf_data[vcf_data["DP"] >= intervals[interval_indx - 1]]

    # 保存过滤后的VCF文件
    dp_filtered_data.to_csv("filtered_by_dp.vcf", sep='\t', index=False)

def filter_by_1kgp3(vcf_data):
    # 读取1KGP3常见变异列表
    common_variants = pd.read_csv("../1kgp3_variants.txt", sep='\t', header=None, names=["CHROM", "POS", "REF", "ALT"])

    # 保留在常见变异列表中的变异
    merged_data = pd.merge(vcf_data, common_variants, on=["CHROM", "POS", "REF", "ALT"])

    # 保存过滤后的VCF文件
    merged_data.to_csv("filtered_by_1kgp3.vcf", sep='\t', index=False)


if __name__ == '__main__' :
    main(sys.argv[1:])

    