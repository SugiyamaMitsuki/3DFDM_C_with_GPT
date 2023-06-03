#ifndef LAYER_PROPERTY_READER_H
#define LAYER_PROPERTY_READER_H

#include <string>
#include <vector>

struct LayerProperty
{
    // 列番号 列名 書式 説明
    // 01 CODE %8c 3次メッシュコード（世界測地系）
    // 02 S0 %.1f 地表標高
    // 03 E0 %.1f 第0層下面標高(m)
    // 04 E1 %.1f 第1層下面標高(m)
    // : : : :
    // 31 E28 %.1f 第28層下面標高(m)
    // 32 E29 %.1f 第29層下面（地震基盤面）標高(m)
    // 33 E30 %.1f 第30層下面標高(m)

    std::string CODE; ///< 3次メッシュコード（世界測地系）
    double S0;
    double E[31];
    // double S0;        ///< 地表標高
    // double E0;        ///< 第0層下面標高(m)
    // double E1;        ///< 第1層下面標高(m)
    // double E2;        ///< 第2層下面標高(m)
    // double E3;        ///< 第3層下面標高(m)
    // double E4;        ///< 第4層下面標高(m)
    // double E5;        ///< 第5層下面標高(m)
    // double E6;        ///< 第6層下面標高(m)
    // double E7;        ///< 第7層下面標高(m)
    // double E8;        ///< 第8層下面標高(m)
    // double E9;        ///< 第9層下面標高(m)
    // double E10;       ///< 第10層下面標高(m)
    // double E11;       ///< 第11層下面標高(m)
    // double E12;       ///< 第12層下面標高(m)
    // double E13;       ///< 第13層下面標高(m)
    // double E14;       ///< 第14層下面標高(m)
    // double E15;       ///< 第15層下面標高(m)
    // double E16;       ///< 第16層下面標高(m)
    // double E17;       ///< 第17層下面標高(m)
    // double E18;       ///< 第18層下面標高(m)
    // double E19;       ///< 第19層下面標高(m)
    // double E20;       ///< 第20層下面標高(m)
    // double E21;       ///< 第21層下面標高(m)
    // double E22;       ///< 第22層下面標高(m)
    // double E23;       ///< 第23層下面標高(m)
    // double E24;       ///< 第24層下面標高(m)
    // double E25;       ///< 第25層下面標高(m)
    // double E26;       ///< 第26層下面標高(m)
    // double E27;       ///< 第27層下面標高(m)
    // double E28;       ///< 第28層下面標高(m)
    // double E29;       ///< 第29層下面（地震基盤面）標高(m)
    // double E30;       ///< 第30層下面標高(m)
    
};

class LayerPropertyReader
{
public:
    LayerPropertyReader(const std::string &filename, double minLat, double maxLat, double minLon, double maxLon); // 範囲指定よみこみ
    LayerProperty operator[](std::string meshcode);                                                               // meshcodeを引数にして、そのメッシュコードのLayer情報を返す
    int getLayerByElevation(LayerProperty targetMesh, double elevation) const; // 標高を引数にして、その標高の物性値番号を返す
    int size() const { return data.size(); }

    std::vector<LayerProperty>  &getData() { return data; }
    
private:
    void loadData();
    LayerProperty parseLine(const std::string& line);
    
    std::string file;
    double minLat, maxLat, minLon, maxLon;
    std::vector<LayerProperty> data;
    LayerProperty findLayerPropertyByCode(const std::string &meshcode);
    LayerProperty findNearestLayerProperty(const std::string &meshcode);

};

#endif // LAYER_PROPERTY_READER_H
