#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <vector>
#include <stdexcept>
#include <sstream>

// 1次元配列クラステンプレート
template <typename T>
class Array1D
{
public:
    // コンストラクタ
    Array1D() : dim{0} {} // デフォルトコンストラクタを追加

    Array1D(int dim, int start_dim = 0) : data(dim), dim{dim}, start{start_dim} {}

    // Constructor for creating Array1D from std::vector
    Array1D(const std::vector<T>& vec, int start_dim = 0)
        : data(vec), dim{static_cast<int>(vec.size())}, start{start_dim} {}
    
    // Constructor for creating Array1D from std::vector
    // Array1D(const std::vector<T>& vec, int start_dim = 0)
    //     : data(vec.begin(), vec.end()), dim{static_cast<int>(vec.size())}, start{start_dim} {}



    // 配列のサイズを設定するメソッド (デフォルトの開始インデックスを使用)
    void resize(int end_dim)
    {
        resize(0, end_dim);
    }

    // 配列のサイズを設定するメソッド (カスタム開始インデックスと終了インデックスを指定)
    void resize(int start_dim, int end_dim)
    {
        int dim = end_dim - start_dim + 1;
        this->dim = dim;
        start = start_dim;
        data.resize(dim);
    }

    // 配列にアクセスするための演算子オーバーロード
    T &operator()(int i)
    {
        return data[index(i)];
    }

    const T &operator()(int i) const
    {
        return data[index(i)];
    }

    // サイズを返す関数
    int size() const
    {
        return dim;
    }

private:
    // 1次元配列へのインデックスを計算する関数
    int index(int i) const
    {
        // if (i < start || i >= start + dim)
        // {
        //     std::stringstream error_message;
        //     error_message << "Array1D index out of range: " << i
        //                   << " is not in [" << start << ", " << start + dim - 1 << "]";
        //     throw std::out_of_range(error_message.str());
        // }
        return i - start;
    }
    std::vector<T> data;
    int dim;
    int start;
};

#endif // ARRAY1D_H
