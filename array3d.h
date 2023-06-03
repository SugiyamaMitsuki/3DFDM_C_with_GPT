#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <vector>
#include <stdexcept>
#include <sstream>

// 3次元配列クラステンプレート
template <typename T>
class Array3D {
public:
    // コンストラクタ
    Array3D() : dims{0, 0, 0} {} // デフォルトコンストラクタを追加

    Array3D(int dim1, int dim2, int dim3, int start_dim1 = 0, int start_dim2 = 0, int start_dim3 = 0) :
        data(dim1 * dim2 * dim3), dims{dim1, dim2, dim3}, start{start_dim1, start_dim2, start_dim3} {}

    // 配列のサイズを設定するメソッド (デフォルトの開始インデックスを使用)
    void resize(int end_dim1, int end_dim2, int end_dim3) {
        resize(0, end_dim1, 0, end_dim2, 0, end_dim3);
    }

    // 配列のサイズを設定するメソッド (カスタム開始インデックスと終了インデックスを指定)
    void resize(int start_dim1, int end_dim1, int start_dim2, int end_dim2, int start_dim3, int end_dim3) {
        int dim1 = end_dim1 - start_dim1 + 1;
        int dim2 = end_dim2 - start_dim2 + 1;
        int dim3 = end_dim3 - start_dim3 + 1;
        dims[0] = dim1;
        dims[1] = dim2;
        dims[2] = dim3;
        start[0] = start_dim1;
        start[1] = start_dim2;
        start[2] = start_dim3;
        data.resize(dim1 * dim2 * dim3);
    }

    // 配列にアクセスするための演算子オーバーロード
    T& operator()(int i, int j, int k) {
        return data[index(i, j, k)];
    }

    const T& operator()(int i, int j, int k) const {
        return data[index(i, j, k)];
    }

    // サイズを返す関数
    int size(int dim) const {
        return dims[dim];
    }

private:
    // 1次元配列へのインデックスを計算する関数
    inline int index(int i, int j, int k) const {
        // if (i < start[0] || i >= start[0] + dims[0] || j < start[1] || j >= start[1] + dims[1] || k < start[2] || k >= start[2] + dims[2]) {
        //     std::stringstream error_message;
        //     error_message << "Array3D index out of range: (" << i << ", " << j << ", " << k << ")"
        //                 << " is not in [" << start[0] << ", " << start[0] + dims[0] - 1 << "] x [" << start[1] << ", " << start[1] + dims[1] - 1 << "] x [" << start[2] << ", " << start[2] + dims[2] - 1 << "]";
        //                 // << ", in function \"" << __FUNCTION__ << "\" at " << __FILE__ << ":" << __LINE__;
        //     throw std::out_of_range(error_message.str());
        // }
        return (i - start[0]) * dims[1] * dims[2] + (j - start[1]) * dims[2] + (k - start[2]);
    }
    std::vector<T> data;
    int dims[3];
    int start[3];

};

#endif // ARRAY3D_H



//使用例

// Array3D<int> a;
//     a.resize(-1, 1, 0, 2, 0, 2); // 3 x 3 x 3 の配列を生成（開始インデックスが -1, 0, 0）


// #include <iostream>
// #include "array3d.h"

// int main() {
//     Array3D<int> a;
//     a.resize(-1, 1, 0, 2, 0, 2); // 3 x 3 x 3 の配列を生成（開始インデックスが -1, 0, 0）

//     // 配列に値を設定
//     for (int i = -1; i <= 1; ++i) {
//         for (int j = 0; j <= 2; ++j) {
//             for (int k = 0; k <= 2; ++k) {
//                 a(i, j, k) = i * 100 + j * 10 + k;
//             }
//         }
//     }

//     // 配列の値を表示
//     for (int i = -1; i <= 1; ++i) {
//         for (int j = 0; j <= 2; ++j) {
//             for (int k = 0; k <= 2; ++k) {
//                 std::cout << "a(" << i << ", " << j << ", " << k << ") = " << a(i, j, k) << std::endl;
//             }
//         }
//     }

//     return 0;
// }
