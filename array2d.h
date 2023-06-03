#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>
#include <stdexcept>
#include <sstream>

// 2次元配列クラステンプレート
template <typename T>
class Array2D
{
public:
    // コンストラクタ
    Array2D() : dims{0, 0} {} // デフォルトコンストラクタを追加

    Array2D(int dim1, int dim2, int start_dim1 = 0, int start_dim2 = 0) : data(dim1 * dim2), dims{dim1, dim2}, start{start_dim1, start_dim2} {}

    // 配列のサイズを設定するメソッド (デフォルトの開始インデックスを使用)
    void resize(int end_dim1, int end_dim2)
    {
        resize(0, end_dim1, 0, end_dim2);
    }

    // 配列のサイズを設定するメソッド (カスタム開始インデックスと終了インデックスを指定)
    void resize(int start_dim1, int end_dim1, int start_dim2, int end_dim2)
    {
        int dim1 = end_dim1 - start_dim1 + 1;
        int dim2 = end_dim2 - start_dim2 + 1;
        dims[0] = dim1;
        dims[1] = dim2;
        start[0] = start_dim1;
        start[1] = start_dim2;
        data.resize(dim1 * dim2);
    }

    // 配列にアクセスするための演算子オーバーロード
    T &operator()(int i, int j)
    {
        return data[index(i, j)];
    }

    const T &operator()(int i, int j) const
    {
        return data[index(i, j)];
    }

    // サイズを返す関数
    int size(int dim) const
    {
        return dims[dim];
    }

private:
    // 1次元配列へのインデックスを計算する関数
    int index(int i, int j) const
    {
        // if (i < start[0] || i >= start[0] + dims[0] || j < start[1] || j >= start[1] + dims[1])
        // {
        //     std::stringstream error_message;
        //     error_message << "Array2D index out of range: (" << i << ", " << j << ")"
        //                   << " is not in [" << start[0] << ", " << start[0] + dims[0] - 1 << "] x [" << start[1] << ", " << start[1] + dims[1] - 1 << "]";
        //     throw std::out_of_range(error_message.str());
        // }
        return (i - start[0]) * dims[1] + (j - start[1]);
    }
    std::vector<T> data;
    int dims[2];
    int start[2];
};

#endif // ARRAY2D_H
