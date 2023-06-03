#ifndef MATERIAL_PROPERTY_READER_H
#define MATERIAL_PROPERTY_READER_H

#include <string>
#include <vector>

struct MaterialProperty {
    int STN;  ///< 物性値番号
    double SVP;  ///< P 波速度(m/s）
    double SVS;  ///< S 波速度(m/s)
    double SRO;  ///< 密度(kg/m3）
    double SQP;  ///< Qp 値※ (1Hz における値)
    double SQS;  ///< Qs 値※ (1Hz における値)
};

class MaterialPropertyReader {
public:
    MaterialPropertyReader(const std::string& filename);
    MaterialProperty operator[](int stn);
    int size() const { return data.size(); }
    

private:
    void loadData();
    std::string file; ///< The name of the CSV file to read.
    std::vector<MaterialProperty> data;
};

#endif // MATERIAL_PROPERTY_READER_H


// https://www.j-shis.bosai.go.jp/map/JSHIS2/data/DOC/DataFileRule/A-RULES.pdf
