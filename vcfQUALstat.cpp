#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>

std::vector<double> get_quals(std::istream &in);
std::vector<int> get_dps(std::istream &in);

template <typename T>
void out2File(const std::string &out, const std::vector<T> &v)
{
    std::ofstream outf(out);
    if (!outf.is_open())
    {
        std::cerr << "Error to open out file." << std::endl;
        return;
    }

    for (const auto &item : v)
    {
        outf << item << std::endl;
    }
    outf.close();
}

int main(int argc, char **argv)
{
    std::ifstream in_vcf_file("data/variants.vcf"); // input vcf file
    if (!in_vcf_file.is_open())
    {
        std::cerr << "Error to open vcf file." << std::endl;
        return 0;
    }

    std::vector<double> QUALs = get_quals(in_vcf_file); // store Qualities of SNP
    in_vcf_file.clear();                                // 清除文件流的状态标志，比如EOF
    in_vcf_file.seekg(0);                               // 将文件指针移动到文件开头
    std::vector<int> DPs = get_dps(in_vcf_file);

    // std::ofstream out("Qualout.txt");

    std::sort(QUALs.begin(), QUALs.end());
    out2File("Qualout.txt", QUALs);
    std::sort(DPs.begin(), DPs.end());
    out2File("DPout.txt", DPs);

    // 计算30%位数
    // size_t percentile_30_index = static_cast<size_t>(0.3 * (QUALs.size() - 1));
    // double percentile_30_value = QUALs[percentile_30_index];

    in_vcf_file.close();
    // out.close();

    return 0;
}

std::vector<double> get_quals(std::istream &in)
{

    std::string line;
    std::vector<double> QUALs;
    double minQualValue = std::numeric_limits<double>::max();
    double maxQualValue = std::numeric_limits<double>::min();
    int count = 0;

    while (std::getline(in, line))
    {
        if (line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string field;
        int column = 0;
        double qualValue;

        while (std::getline(ss, field, '\t'))
        {
            column++;
            if (column == 6)
            { // QUAL 在第6列
                try
                {
                    qualValue = std::stod(field);
                }
                catch (const std::invalid_argument &e)
                {
                    std::cerr << "Invalid qual value" << field << std::endl;
                    continue;
                }

                QUALs.push_back(qualValue);
                if (qualValue > maxQualValue)
                    maxQualValue = qualValue;
                if (qualValue < minQualValue)
                    minQualValue = qualValue;
                count++;
            }
        }
    }

    return QUALs;
}

std::vector<int> get_dps(std::istream &in)
{
    std::string line;
    std::vector<int> DPs;

    while (std::getline(in, line))
    {
        if (line[0] == '#')
            continue;

        std::stringstream ss(line);
        std::string fields;
        int column = 0;

        while (std::getline(ss, fields, '\t'))
        {
            column++;
            if (column == 8)
            {
                std::stringstream INFOs(fields);
                std::string info;
                int dpValue;

                while (std::getline(INFOs, info, ';'))
                {
                    if (info.substr(0, 3) == "DP=") // 检查是否是 DP 字段
                    {
                        try
                        {
                            int dpValue = std::stoi(info.substr(3)); // 跳过 "DP=" 提取值
                            // std::cout << dpValue << std::endl;
                            DPs.push_back(dpValue);
                        }
                        catch (const std::invalid_argument &e)
                        {
                            std::cerr << "Invalid DP value: " << info << std::endl;
                            continue;
                        }
                    }
                }
            }
        }
    }
    return DPs;
}
