#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <numeric> 
#include <chrono>
#include <functional>
#include <cmath>
#include <queue>
#include <unordered_map>
using namespace std;

template <class T>
class SparseMatrix; //前向声明

template <class T>
class Matrix {
private:
    vector<vector<T>> data;
    size_t rows=0;
    size_t cols=0;

public:
    //初始化
    //默认构造函数
    Matrix() = default;
    //显式构造函数
    Matrix(size_t r, size_t c) : rows(r), cols(c), data(r, vector<T>(c)) {}

    //向量初始化
    Matrix(const vector<T>& rowVector) : rows(1), cols(rowVector.size()), data(1, rowVector) {}

    SparseMatrix<T> toSparse() const;

    //列表初始化法
    Matrix(initializer_list<initializer_list<T>> init) {
        rows = init.size();
        cols = init.begin()->size();
        data.resize(rows);
        auto data_ptr = data.begin();
        for (auto row = init.begin(); row != init.end(); ++row, ++data_ptr) {
            data_ptr->assign(row->begin(), row->end());
        }
    }

    //矩阵乘积
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            cout << cols << "!=" << other.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;//精度控制
            }
        }
        return result;
    }

    //if优化乘法
    static Matrix<T> if_mul(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < a.cols; ++k) {
                    if (a.data[i][k] == 0 || b.data[k][j] == 0) result.data[i][j] = 0;
                    else result.data[i][j] += a.data[i][k] * b.data[k][j];
                }
                // 精度控制
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;
            }
        }
        return result;
    }

    //索引优化乘法(左乘)
    static Matrix<T> ind_mul_l(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i+=2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //对块内每行
                    for (size_t l = 0; l < 2; ++l) {  //对块内每列
                        result.data[i + k][j] += a.data[i + k][i + l] * b.data[i + l][j];
                    }
                    //if (abs(result.data[i + k][j]) < 1e-11) result.data[i + k][j] = 0;  //精度控制
                }
            }
        }
        return result;
    }

    //索引优化乘法(右乘)
    static Matrix<T> ind_mul_r(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < result.cols; j += 2) {
                for (size_t k = 0; k < a.cols; ++k) {
                    for (size_t l = 0; l < 2; ++l) {
                        result.data[i][j + l] += a.data[i][k] * b.data[k][j + l];
                    }
                }
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;  //精度控制
            }
        }
        return result;
    }

    //索引优化乘法(左乘，捎带置换，先乘后置换)
    static Matrix<T> ind_mul_left(const Matrix& a, const Matrix& b, vector<int> p) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i += 2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //对块内每行
                    for (size_t l = 0; l < 2; ++l) {  //对块内每列
                        result.data[p[i + k]][j] += a.data[i + k][i + l] * b.data[i + l][j];
                    }
                }
                //if (abs(result.data[p[i]][j]) < 1e-11) result.data[p[i]][j] = 0;  //精度控制
            }
        }
        return result;
    }

    //索引优化乘法(左乘，捎带置换，先置换后乘)
    static Matrix<T> ind_mul_left_recover(const Matrix& a, const Matrix& b, vector<int> p) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i += 2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //对块内每行
                    for (size_t l = 0; l < 2; ++l) {  //对块内每列
                        result.data[i + k][j] += a.data[i + k][i + l] * b.data[p[i + l]][j];
                    }
                }
                //if (abs(result.data[p[i]][j]) < 1e-11) result.data[p[i]][j] = 0;  //精度控制
            }
        }
        return result;
    }

    //索引优化乘法(右乘，捎带置换，先乘后置换)
    static Matrix<T> ind_mul_right(const Matrix& a, const Matrix& b, const vector<int>& p) {
        // 检查矩阵乘法的维度是否匹配
        if (a.cols != b.rows) {
            cout << a.cols << " != " << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols); 
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < b.cols; j += 2) {
                for (size_t k = 0; k < 2 ; ++k) { // 对块内每列
                    for (size_t l = 0; l < 2; ++l) { // 对块内每行
                        result.data[i][p[j + k]] += a.data[i][j + l] * b.data[j + l][j + k];
                    }
                    // 精度控制
                    // if (abs(result.data[i][p[j + k]]) < 1e-11) result.data[i][p[j + k]] = 0;
                }
            }
        }

        return result;
    }

    //索引优化乘法(右乘，捎带置换，先置换后乘)
    static Matrix<T> ind_mul_right_recover(const Matrix& a, const Matrix& b, const vector<int>& p) {
        // 检查矩阵乘法的维度是否匹配
        if (a.cols != b.rows) {
            cout << a.cols << " != " << b.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < b.cols; j += 2) {
                for (size_t k = 0; k < 2; ++k) { // 对块内每列
                    for (size_t l = 0; l < 2; ++l) { // 对块内每行
                        result.data[i][j + k] += a.data[i][p[j + l]] * b.data[j + l][j + k];
                    }
                    // 精度控制
                    // if (abs(result.data[i][p[j + k]]) < 1e-11) result.data[i][p[j + k]] = 0;
                }
            }
        }

        return result;
    }
  

    Matrix<T> mulByr(const Matrix<T>& other) const {
        if (cols != other.getRows() || other.getCols() != 1) {
            cout << cols << "!=" << other.rows << endl;
            throw invalid_argument("乘法列与行不匹配");
        }

        Matrix<T> result(rows, 1);
        for (size_t j = 0; j < cols; ++j) {
            if (other.data[j][0] == 1) {
                for (size_t i = 0; i < rows; ++i) {
                    result.data[i][0] += data[i][j];  //累加A的第j列到结果矩阵
                }
            }
        }
        return result;
    }

    //矩阵数乘
    Matrix operator*(const T& scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    //矩阵加法
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("加法行列不匹配");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    //矩阵减法
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("减法行列不匹配");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    //矩阵转置
    Matrix<T> transpose() const {
        Matrix<T> result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    //矩阵的逆
    Matrix<T> inverse() const {
        if (rows != cols) {
            throw invalid_argument("矩阵必须是方阵");
        }
        Matrix<T> result(rows, cols);

        if (rows == 2) {
            result[0][0] = data[1][1];
            result[1][1] = data[0][0];
            result[0][1] = -data[0][1];
            result[1][0] = -data[1][0];
            return result * (1 / (data[0][0] * data[1][1] - data[0][1] * data[1][0]));
        }

        Matrix<T> copy(*this);

        for (size_t i = 0; i < rows; ++i) { //单位矩阵
            result[i][i] = static_cast<T>(1);
        }

        for (size_t i = 0; i < rows; ++i) {
            // 寻找主元
            if (copy[i][i] == static_cast<T>(0)) {
                bool found = false;
                for (size_t j = i + 1; j < rows && !found; ++j) {
                    if (copy[j][i] != static_cast<T>(0)) {
                        // 交换行
                        swap(copy.data[i], copy.data[j]);
                        swap(result.data[i], result.data[j]);
                        found = true;
                    }
                }
                if (!found) {
                    throw runtime_error("矩阵不可逆");
                }
            }

            // 使主元变为1
            T div = copy[i][i];
            for (size_t j = 0; j < cols; ++j) {
                copy[i][j] /= div;
                result[i][j] /= div;
            }

            // 消去其他行的该列元素
            for (size_t k = 0; k < rows; ++k) {
                if (k != i && copy[k][i] != static_cast<T>(0)) {
                    T factor = copy[k][i];
                    for (size_t j = 0; j < cols; ++j) {
                        copy[k][j] -= copy[i][j] * factor;
                        result[k][j] -= result[i][j] * factor;
                    }
                }
            }
        }

        return result;
    }

    //矩阵访问
    vector<T>& operator[](size_t index) {
        if (index >= rows) {
            throw std::out_of_range("访问越界");
        }
        return data[index];
    }

    //矩阵等值比较
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) return false;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (abs(data[i][j] - other.data[i][j]) > 1e-3) return false;
            }
        }
        return true;
    }

    //矩阵第二范式单位化(||A||)
    static void normalizeRows(Matrix<T>& matrix) {
        for (size_t i = 0; i < matrix.getRows(); ++i) {
            T norm = 0;
            for (size_t j = 0; j < matrix.getCols(); ++j) {
                norm += matrix[i][j] * matrix[i][j];
            }
            norm = sqrt(norm);
            if (norm > 0) { // 防止除以零
                for (size_t j = 0; j < matrix.getCols(); ++j) {
                    matrix[i][j] /= norm;
                }
            }
        }
    }

    //打印
    void print() const {
        for (const auto& row : data) {
            for (const auto& elem : row) {
                cout << elem << " ";
            }
            cout << endl;
        }
    }

    //获取行数
    size_t getRows() const { return rows; }

    //获取列数
    size_t getCols() const { return cols; }

    //获取单行作为向量
    Matrix<T> getRow(size_t index) const {
        if (index >= rows) {
            throw std::out_of_range("访问越界");
        }
        return Matrix<T>{data[index]};
    }

    //欧氏距离
    static T distance(Matrix a, Matrix b) {
        if (a.getCols() != b.getCols() || a.getRows() != 1 || b.getRows() != 1) {
            throw std::invalid_argument("行列数不匹配");
        }
        T dist = 0;
        for (size_t i = 0; i < a.getCols(); i++) {
            T diff = a[0][i] - b[0][i];
            dist += diff * diff;
        }
        return sqrt(dist);
    }

    //余弦相似度
    static Matrix<T> cos(Matrix<T> A, Matrix<T> B) {
        Matrix Anorm = A, Bnorm = B;
        //单位化
        normalizeRows(Anorm);
        normalizeRows(Bnorm);

        return Anorm * Bnorm.transpose();
    }

    //生成随机矩阵
    static Matrix random(size_t r, size_t c, T min = numeric_limits<T>::min(), T max = numeric_limits<T>::max()) {
        random_device rd;
        mt19937 gen(rd()); //伪随机数生成器
        uniform_real_distribution<> dis(min, max);

        Matrix randomMatrix(r, c);
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < c; ++j) {
                //如果是整数类型则需要转换
                if constexpr (std::is_integral_v<T>) {
                    randomMatrix[i][j] = static_cast<T>(dis(gen));
                }
                else {
                    randomMatrix[i][j] = dis(gen);
                }
            }
        }
        return randomMatrix;
    }

    //Householder变换
    static Matrix<T> householder(size_t n, T min, T max) {
        Matrix<T> result(n, n);
        if (n % 2 != 0)  throw runtime_error("矩阵不是2的倍数");
        for (int i = 0; i < n / 2; i++) {
            Matrix<T> v = random(2, 1, 0, 10);
            Matrix<T> I = {
                {1,0},
                {0,1}
            };
            T a = (v.transpose() * v)[0][0];
            Matrix<T> H = I - v * v.transpose() * (2 / a);

            result[i * 2][i * 2] = H[0][0];
            result[i * 2][i * 2 + 1] = H[0][1];
            result[i * 2 + 1][i * 2] = H[1][0];
            result[i * 2 + 1][i * 2 + 1] = H[1][1];
        }
        return result;
    }

    //随机生成置换矩阵
    static Matrix<T> RPMatrix(size_t size) {
        Matrix<T> PMatrix(size, size);

        //填充单位矩阵
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[i][i] = static_cast<T>(1);
        }

        //随机数生成器
        random_device rd;
        mt19937 generator(rd());

        vector<int> indices(size);
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), generator);
        //置换
        for (size_t i = 0; i < size; ++i) {
            if (i != indices[i]) {
                swap(PMatrix.data[i], PMatrix.data[indices[i]]);
            }
        }
        return PMatrix;
    }

    //随机置换映射
    static vector<int> RPermutation(size_t size) {
        vector<int> i(size);
        iota(i.begin(), i.end(), 0);
        random_device rd;
        mt19937 generator(rd());
        shuffle(i.begin(), i.end(), generator);

        return i;
    }

    //随机置换矩阵，左乘，交换行(P*A)
    static Matrix<T> leftPMatrix(const vector<int>& p) {
        size_t size = p.size();
        Matrix<T> PMatrix(size, size);
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[p[i]][i] = static_cast<T>(1);  //将1放在置换指定行
        }
        return PMatrix;
    }

    //随机置换矩阵，右乘，交换列(A*P)
    static Matrix<T> rightPMatrix(const vector<int>& p) {
        size_t size = p.size();
        Matrix<T> PMatrix(size, size);
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[i][p[i]] = static_cast<T>(1);  //将1放在置换指定列
        }
        return PMatrix;
    }

    //生成随机验证向量
    static Matrix<T> R01r(size_t length) {
        Matrix<T> vec(length, 1);
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> distrib(0, 1);
        for (size_t i = 0; i < length; ++i) {
            vec.data[i][0] = distrib(gen);
        }
        return vec;
    }


    //strassen算法
    static Matrix<T> Strassen(Matrix<T> A, Matrix<T> B, int &ss, int &cc) {
        if (A.getCols() % 2 != 0 || B.getCols() % 2 != 0) throw runtime_error("矩阵不是2的倍数");
        T n = A.getRows() / 2;
        T m = A.getCols() / 2;
        T k = B.getCols() / 2;


        Matrix<T> A11(n, m);
        Matrix<T> A12(n, m);
        Matrix<T> A21(n, m);
        Matrix<T> A22(n, m);
        Matrix<T> B11(m, k);
        Matrix<T> B12(m, k);
        Matrix<T> B21(m, k);
        Matrix<T> B22(m, k);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A11[i][j] = A[i][j];
                A21[i][j] = A[i + n][j];
                A12[i][j] = A[i][j + m];
                A22[i][j] = A[i + n][j + m];
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                B11[i][j] = B[i][j];
                B21[i][j] = B[i + m][j];
                B12[i][j] = B[i][j + k];
                B22[i][j] = B[i + m][j + k];
            }
        }

        auto start1 = chrono::high_resolution_clock::now();
        Matrix<T> R1 = B12 - B22;
        Matrix<T> R2 = A11 + A12;
        Matrix<T> R3 = A21 + A22;
        Matrix<T> R4 = B21 - B11;
        Matrix<T> R5 = A11 + A22;
        Matrix<T> R6 = B11 - B22;
        Matrix<T> R7 = A12 - A22;
        Matrix<T> R8 = B21 + B22;
        Matrix<T> R9 = A11 - A21;
        Matrix<T> R10 = B11 - B12;
        auto end1 = chrono::high_resolution_clock::now();
        auto t1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1);
        // << "客户端加法预处理计算用时:" << " " << t1.count() << endl;
        ss += t1.count();
        //服务器完成的部分----------------------------------------
        auto start2 = chrono::high_resolution_clock::now();
        Matrix<T> M1 = A11 * R1;
        Matrix<T> M2 = R2 * B22;
        auto end2 = chrono::high_resolution_clock::now();
        auto t2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        //cout << "S1计算用时:" << " " << t2.count() << endl;

        auto start3 = chrono::high_resolution_clock::now();
        Matrix<T> M3 = R3 * B11;
        Matrix<T> M4 = A22 * R4;
        auto end3 = chrono::high_resolution_clock::now();
        auto t3 = chrono::duration_cast<chrono::milliseconds>(end3 - start3);
        //cout << "S2计算用时:" << " " << t3.count() << endl;

        auto start4 = chrono::high_resolution_clock::now();
        Matrix<T> M5 = R5 * R6;
        Matrix<T> M6 = R7 * R8;
        auto end4 = chrono::high_resolution_clock::now();
        auto t4 = chrono::duration_cast<chrono::milliseconds>(end4 - start4);
        //cout << "S3计算用时:" << " " << t4.count() << endl;

        auto start5 = chrono::high_resolution_clock::now();
        Matrix<T> M7 = R9* R10;    
        auto end5 = chrono::high_resolution_clock::now();
        auto t5 = chrono::duration_cast<chrono::milliseconds>(end5 - start5);
        //cout << "S4计算用时:" << " " << t5.count() << endl;

        auto tt = t2 + t3 + t4 + t5 ;
        //cout << "strassen方法云服务器总用时:" << " " << tt.count() << endl;
        ss += tt.count();

        //--------------------------------------------------------
        auto start9 = chrono::high_resolution_clock::now();
        Matrix<T> U11 = M5 + M4 - M2 + M6;
        Matrix<T> U12 = M1 + M2;
        Matrix<T> U21 = M3 + M4;
        Matrix<T> U22 = M5 + M1 - M3 - M7;
        auto end9 = chrono::high_resolution_clock::now();
        auto t9 = chrono::duration_cast<chrono::milliseconds>(end9 - start9);
        //cout << "客户端加法后处理计算用时:" << " " << t9.count() << endl;
        cc += t9.count();
        cout << endl;

        Matrix<T> C(2 * n, 2 * k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                C[i][j] = U11[i][j];
                C[i + n][j] = U12[i][j];
                C[i][j + k] = U21[i][j];
                C[i + n][j + k] = U22[i][j];
            }
        }
        return C;
    }

    //Coppersmith-Winograd算法(一轮，客户端承担更多加法)
    static Matrix<T> Coppersmith(Matrix<T> A, Matrix<T> B, int &s, int &c) {
        if (A.getCols() % 2 != 0 || B.getCols() % 2 != 0) throw runtime_error("矩阵不是2的倍数");
        T n = A.getRows() / 2;
        T m = A.getCols() / 2;
        T k = B.getCols() / 2;


        Matrix<T> A11(n, m);
        Matrix<T> A12(n, m);
        Matrix<T> A21(n, m);
        Matrix<T> A22(n, m);
        Matrix<T> B11(m, k);
        Matrix<T> B12(m, k);
        Matrix<T> B21(m, k);
        Matrix<T> B22(m, k);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A11[i][j] = A[i][j];
                A21[i][j] = A[i + n][j];
                A12[i][j] = A[i][j + m];
                A22[i][j] = A[i + n][j + m];
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < k; j++) {
                B11[i][j] = B[i][j];
                B21[i][j] = B[i + m][j];
                B12[i][j] = B[i][j + k];
                B22[i][j] = B[i + m][j + k];
            }
        }

        auto start1 = chrono::high_resolution_clock::now();
        Matrix<T> S1 = A21 + A22;
        Matrix<T> S2 = S1 - A11;
        Matrix<T> T1 = B12 - B11;
        Matrix<T> T2 = B22 - T1;
        auto end1 = chrono::high_resolution_clock::now();
        auto t1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1);
        // << "客户端加法预处理计算用时:" << " " << t1.count() << endl;
        c += t1.count();
        //服务器完成的部分----------------------------------------
        auto start2 = chrono::high_resolution_clock::now();
        Matrix<T> M1 = A11 * B11;
        auto end2 = chrono::high_resolution_clock::now();
        auto t2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        //cout << "S1计算用时:" << " " << t2.count() << endl;

        auto start3 = chrono::high_resolution_clock::now();
        Matrix<T> M2 = A12 * B21;
        auto end3 = chrono::high_resolution_clock::now();
        auto t3 = chrono::duration_cast<chrono::milliseconds>(end3 - start3);
        //cout << "S2计算用时:" << " " << t3.count() << endl;

        auto start4 = chrono::high_resolution_clock::now();
        Matrix<T> S4 = A12 - S2;
        Matrix<T> M3 = S4 * B22;
        auto end4 = chrono::high_resolution_clock::now();
        auto t4 = chrono::duration_cast<chrono::milliseconds>(end4 - start4);
        //cout << "S3计算用时:" << " " << t4.count() << endl;

        auto start5 = chrono::high_resolution_clock::now();
        Matrix<T> T4 = T2 - B21;
        Matrix<T> M4 = A22 * T4;
        auto end5 = chrono::high_resolution_clock::now();
        auto t5 = chrono::duration_cast<chrono::milliseconds>(end5 - start5);
        //cout << "S4计算用时:" << " " << t5.count() << endl;

        auto start6 = chrono::high_resolution_clock::now();
        Matrix<T> M5 = S1 * T1;
        auto end6 = chrono::high_resolution_clock::now();
        auto t6 = chrono::duration_cast<chrono::milliseconds>(end6 - start6);
        //cout << "S5计算用时:" << " " << t6.count() << endl;

        auto start7 = chrono::high_resolution_clock::now();
        Matrix<T> M6 = S2 * T2;
        auto end7 = chrono::high_resolution_clock::now();
        auto t7 = chrono::duration_cast<chrono::milliseconds>(end7 - start7);
        //cout << "S6计算用时:" << " " << t7.count() << endl;

        auto start8 = chrono::high_resolution_clock::now();
        Matrix<T> S3 = A11 - A21;
        Matrix<T> T3 = B22 - B12;
        Matrix<T> M7 = S3 * T3;
        auto end8 = chrono::high_resolution_clock::now();
        auto t8 = chrono::duration_cast<chrono::milliseconds>(end8 - start8);
        //cout << "S7计算用时:" << " " << t8.count() << endl;

        auto tt = t2 + t3 + t4 + t5 + t6 + t7 + t8;
        //cout << "云服务器总用时:" << " " << tt.count() << endl;
        s += tt.count();

        //--------------------------------------------------------
        auto start9 = chrono::high_resolution_clock::now();
        Matrix<T> U1 = M1 + M2;
        Matrix<T> U2 = M1 + M6;
        Matrix<T> U3 = U2 + M7;
        Matrix<T> U4 = U2 + M5;
        Matrix<T> U5 = U4 + M3;
        Matrix<T> U6 = U3 - M4;
        Matrix<T> U7 = U3 + M5;
        auto end9 = chrono::high_resolution_clock::now();
        auto t9 = chrono::duration_cast<chrono::milliseconds>(end9 - start9);
        //cout << "客户端加法后处理计算用时:" << " " << t9.count() << endl;
        c += t9.count();
        cout << endl;

        Matrix<T> C(2 * n, 2 * k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                C[i][j] = U1[i][j];
                C[i + n][j] = U6[i][j];
                C[i][j + k] = U5[i][j];
                C[i + n][j + k] = U7[i][j];
            }
        }
        return C;
    }
};

#include "SparseMatrix.h"

//将普通矩阵转换为稀疏矩阵
template <class T>
SparseMatrix<T> Matrix<T>::toSparse() const {
    SparseMatrix<T> sparse(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (data[i][j] != 0) {
                sparse.set(i, j, data[i][j]);
            }
        }
    }
    return sparse;
}















//K-means聚类安全外包
template <typename T>
pair<Matrix<T>, Matrix<T>> kMeans(Matrix<T>& points, size_t k, int maxIterations) {
    chrono::high_resolution_clock::time_point start;
    start = chrono::high_resolution_clock::now();
    chrono::high_resolution_clock::time_point end;
    end = chrono::high_resolution_clock::now();
    chrono::milliseconds t;
    chrono::milliseconds blitime;
    chrono::milliseconds vertime;
    chrono::milliseconds rectime;
    int liutime;
    int s = 0;
    int c = 0;
    int ss = 0;
    int cc = 0;
    blitime = chrono::duration_cast<chrono::milliseconds>(start - start);
    vertime = chrono::duration_cast<chrono::milliseconds>(start - start);
    rectime = chrono::duration_cast<chrono::milliseconds>(start - start);
    
    size_t v = points.getRows(), d = points.getCols();
    Matrix<T> centroids(k, d); //质心矩阵
    vector<int> labels(v, 0);
    int sp = 10;

    //第一步：初始化，随机选择质心
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, v - 1);
    for (size_t i = 0; i < k; ++i) {
        size_t randIdx = dis(gen);
        for (size_t j = 0; j < d; ++j) {
            centroids[i][j] = points[randIdx][j];
        }
    }

    //第二步：单位化
    Matrix<T>::normalizeRows(points);
    Matrix<T>::normalizeRows(centroids);
    Matrix<T> W_D = points;
    Matrix<T> W_T = centroids.transpose();

    int iter = 0;
    for (; iter < 1; ++iter) {
        bool changed = false;

        //第三步：密钥生成
        Matrix<double> Q1 = Matrix<double>::householder(v, -pow(2, sp), pow(2, sp));
        Matrix<double> Q2 = Matrix<double>::householder(d, -pow(2, sp), pow(2, sp));
        Matrix<double> Q3 = Matrix<double>::householder(k, -pow(2, sp), pow(2, sp));
        auto p1 = Matrix<double>::RPermutation(v);
        auto p2 = Matrix<double>::RPermutation(d);
        auto p3 = Matrix<double>::RPermutation(k);

        //第四步：盲化
        start = chrono::high_resolution_clock::now();
        Matrix<T> W_D_t = Matrix<double>::ind_mul_left(Q1, W_D, p1);
        Matrix<T> W_T_t = Matrix<double>::ind_mul_left(Q2, W_T, p2);
        Matrix<T> W_D_b = Matrix<double>::ind_mul_right(W_D_t, Q2, p2);
        Matrix<T> W_T_b = Matrix<double>::ind_mul_right(W_T_t, Q3, p3);
        Matrix<T> r[10];
        Matrix<double> V[10];
        for (int i = 0; i < 10; i++) {
            r[i] = Matrix<T>::R01r(k);
            V[i] = W_T_b.mulByr(r[i]);//用于验证
        }
        end = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end - start);
        blitime += t;
        cout << "盲化用时:" << " " << t.count() << endl;

        //第五步：计算
        Matrix<T> C_b = Matrix<double>::Coppersmith(W_D_b, W_T_b, s, c);

        Matrix<T> t1 = Matrix<double>::Strassen(W_D_b, W_T_b, ss, cc);

        //第六步：验证  
        t = chrono::duration_cast<chrono::milliseconds>(start - start);
        Matrix<double> temp,temp2;
        for (int i = 0; i < 10; i++) {
            start = chrono::high_resolution_clock::now();
            temp = C_b.mulByr(r[i]);
            temp2 = W_D_b * V[i];
            end = chrono::high_resolution_clock::now();
            t += chrono::duration_cast<chrono::milliseconds>(end - start);
            if (!(temp == temp2)) {
                cout << "第" << i + 1 << "轮验证未通过" << endl;
                break;
            }
            else cout << "第" << i + 1 << "轮验证通过" << endl;

        }
        vertime += t;
        cout << "验证用时:" << " " << t.count() << endl;
       
        //第七步：恢复
        start = chrono::high_resolution_clock::now();
        Matrix<T> C_h = Matrix<double>::ind_mul_left_recover(Q1, C_b, p1);
        Matrix<T> C = Matrix<double>::ind_mul_right_recover(C_h, Q3, p3);
        end = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end - start);
        rectime += t;
        cout << "恢复用时:" << " " << t.count() << endl;
        cout << endl;

        //第八步：更新
        Matrix<T> similarity = C;
        // 重新分配最相似的质心
        for (size_t i = 0; i < v; ++i) {
            T maxSimilarity = -numeric_limits<T>::max();
            int bestCluster = 0;
            for (size_t j = 0; j < k; ++j) {
                if (similarity[i][j] > maxSimilarity) {
                    maxSimilarity = similarity[i][j];
                    bestCluster = j;
                }
            }
            if (labels[i] != bestCluster) {
                labels[i] = bestCluster;
                changed = true;
            }
        }

        // 更新质心
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < d; j++) {
                centroids[i][j] = 0;
            }
        }
        vector<int> count(k, 0);
        for (size_t i = 0; i < v; ++i) {
            count[labels[i]]++;
            for (size_t j = 0; j < d; ++j) {
                centroids[labels[i]][j] += points[i][j];
            }
        }
        for (size_t j = 0; j < k; ++j) {
            if (count[j] > 0) {
                for (size_t l = 0; l < d; ++l) {
                    centroids[j][l] /= count[j];
                }
            }
        }
        Matrix<T>::normalizeRows(centroids);
        W_T = centroids.transpose();
        //if (!changed) break; //如果没有变化，提前终止


        Matrix<T> r1 = Matrix<T>::R01r(d);
        Matrix<T> r2 = Matrix<T>::R01r(k);
        start = chrono::high_resolution_clock::now();
        Matrix<T> tt1 = W_D_t * r1;
        Matrix<T> tt2 = W_T_t * r2;
        Matrix<T> tt3 = W_D_t + W_D_t + W_D_t;
        Matrix<T> tt4 = W_T_t + W_T_t + W_T_t;
        Matrix<T> tt5 = C_h * r2;
        Matrix<T> tt6 = C_h + C_h + C_h;
        end = chrono::high_resolution_clock::now();
        t= chrono::duration_cast<chrono::milliseconds>(end - start);
        liutime = blitime.count() + vertime.count() + rectime.count() + t.count();
    }

    cout << "迭代次数： " << iter*5  << endl;
    cout << endl;
    /*
    vector<queue<int>> queues(k);
    for (size_t i = 0; i < v; ++i) {
        queues[labels[i]].push(i);
    }
    for (size_t i = 0; i < k; i++) {
        cout << "簇" << i << "包含数据点：";
        while (!queues[i].empty()) {
            cout << queues[i].front() << " ";
            queues[i].pop();
        }
        cout << endl;
    }
    */
    cout << endl;
    cout << "盲化k-means总用时：" << blitime.count()*5 << endl;
    cout << "验证k-means总用时：" << vertime.count()*5 << endl;
    cout << "恢复k-means总用时：" << rectime.count()*5 << endl;
    cout << "云服务器k-means总用时：" << s * 5 << endl;
    cout << "客户端k-means预处理用时：" << c * 5 << endl;
    cout << endl;
    cout << "strassen的k-means总用时：" << ss * 5 << endl;
    cout << endl;
    cout << "liu方法k-means总用时：" << liutime * 5 << endl;
    cout << endl;

    return { points, centroids };
}


//连乘安全外包
template <typename T>
Matrix<T> cMultiply(Matrix<T>& X, Matrix<T>& A, Matrix<T>& B) {
    chrono::high_resolution_clock::time_point start;
    start = chrono::high_resolution_clock::now();
    chrono::high_resolution_clock::time_point end;
    end = chrono::high_resolution_clock::now();
    chrono::milliseconds t;
    int s = 0;
    int c = 0;
    int ss = 0;
    int cc = 0;
    int liutime = 0;

    size_t n = X.getRows(), v = A.getRows(), d = A.getCols(), k = B.getRows();

    int sp = 10;

    //第一步：单位化
    Matrix<T>::normalizeRows(A);
    Matrix<T>::normalizeRows(B);
    Matrix<T> W_D = A;
    Matrix<T> W_T = B.transpose();

    //第二步：密钥生成
    Matrix<double> Q1 = Matrix<double>::householder(n, -pow(2, sp), pow(2, sp));
    Matrix<double> Q2 = Matrix<double>::householder(v, -pow(2, sp), pow(2, sp));
    Matrix<double> Q3 = Matrix<double>::householder(d, -pow(2, sp), pow(2, sp));
    Matrix<double> Q4 = Matrix<double>::householder(k, -pow(2, sp), pow(2, sp));
    auto p1 = Matrix<double>::RPermutation(n);
    auto p2 = Matrix<double>::RPermutation(v);
    auto p3 = Matrix<double>::RPermutation(d);
    auto p4 = Matrix<double>::RPermutation(k);

    //第三步：盲化
    start = chrono::high_resolution_clock::now();
    Matrix<T> X_t = Matrix<double>::ind_mul_left(Q1, X, p1);
    Matrix<T> W_D_t = Matrix<double>::ind_mul_left(Q2, W_D, p2);
    Matrix<T> W_T_t = Matrix<double>::ind_mul_left(Q3, W_T, p3);
    Matrix<T> X_b = Matrix<double>::ind_mul_right(X_t, Q2, p2);
    Matrix<T> W_D_b = Matrix<double>::ind_mul_right(W_D_t, Q3, p3);
    Matrix<T> W_T_b = Matrix<double>::ind_mul_right(W_T_t, Q4, p4);
    Matrix<T> r[10];
    Matrix<double> V[10];
    for (int i = 0; i < 10; i++) {
        r[i] = Matrix<T>::R01r(k);
        V[i] = W_T_b.mulByr(r[i]);//用于验证
    }
    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    liutime += t.count();
    cout << "盲化用时:" << " " << t.count() << endl;

    //第四步：计算
    Matrix<T> M = Matrix<double>::Coppersmith(W_D_b, W_T_b, s ,c);
    Matrix<T> Z_b = Matrix<double>::Coppersmith(X_b, M, s, c);

    Matrix<T> t1 = Matrix<double>::Strassen(W_D_b, W_T_b, ss, cc);
    Matrix<T> t2 = Matrix<double>::Strassen(X_b, t1, ss, cc);

    //第五步：验证
    t = chrono::duration_cast<chrono::milliseconds>(start - start);
    Matrix<double> temp, temp2;
    for (int i = 0; i < 10; i++) {
        start = chrono::high_resolution_clock::now();
        temp = Z_b.mulByr(r[i]);
        temp2 = X_b * (W_D_b * V[i]);
        end = chrono::high_resolution_clock::now();
        t += chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!(temp == temp2)) {
            cout << "第" << i + 1 << "轮验证未通过" << endl;
            break;
        }
        else cout << "第" << i + 1 << "轮验证通过" << endl;

    }
    liutime += t.count();
    cout << "验证用时:" << " " << t.count() << endl;

    //第六步：恢复
    start = chrono::high_resolution_clock::now();
    Matrix<T> Z_h = Matrix<double>::ind_mul_left_recover(Q1, Z_b, p1);
    Matrix<T> Z = Matrix<double>::ind_mul_right_recover(Z_h, Q4, p4);
    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    liutime += t.count();
    cout << "恢复用时:" << " " << t.count() << endl;
    cout << endl;


    Matrix<T> r1 = Matrix<T>::R01r(v);
    Matrix<T> r2 = Matrix<T>::R01r(d);
    Matrix<T> r3 = Matrix<T>::R01r(k);
    start = chrono::high_resolution_clock::now();
    Matrix<T> tt1 = X_t * r1;
    Matrix<T> tt2 = W_D_t * r2;
    Matrix<T> tt3 = W_T_t * r3;
    Matrix<T> tt4 = X_t + X_t + X_t;
    Matrix<T> tt5 = W_D_t + W_D_t + W_D_t;
    Matrix<T> tt6 = W_T_t + W_T_t + W_T_t;
    Matrix<T> tt7 = Z_h * r3;
    Matrix<T> tt8 = Z_h + Z_h + Z_h;
    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    liutime += t.count();


    cout << "云服务器连乘用时：" << s << endl;
    cout << "客户端预处理用时：" << c << endl;
    cout << endl;
    cout << "strassen连乘用时：" << ss << endl;
    cout << endl;
    cout << "liu方法连乘总用时：" << liutime << endl;
    cout << endl;

    return Z;
}






//k-means原始计算
template <typename T>
pair<Matrix<T>, Matrix<T>> kMeans_ori(Matrix<T>& points, size_t k, int maxIterations) {
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    chrono::high_resolution_clock::time_point start2;
    chrono::high_resolution_clock::time_point end2;
    chrono::milliseconds t;

    size_t v = points.getRows(), d = points.getCols();
    Matrix<T> centroids(k, d); //质心矩阵
    vector<int> labels(v, 0);

    start = chrono::high_resolution_clock::now();
    //第一步：初始化，随机选择质心
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, v - 1);
    for (size_t i = 0; i < k; ++i) {
        size_t randIdx = dis(gen);
        for (size_t j = 0; j < d; ++j) {
            centroids[i][j] = points[randIdx][j];
        }
    }

    //第二步：单位化
    Matrix<T>::normalizeRows(points);
    Matrix<T>::normalizeRows(centroids);
    Matrix<T> W_D = points;
    Matrix<T> W_T = centroids.transpose();

    int iter = 0;
    for (; iter < 1; ++iter) {
        bool changed = false;

        //第三步：计算
        start2 = chrono::high_resolution_clock::now();
        Matrix<T> C = W_D * W_T;
        end2 = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        cout << "k-means直接计算用时:" << " " << t.count() * 5 << endl;
        cout << endl;

        //第四步：更新
        Matrix<T> similarity = C;
        // 重新分配最相似的质心
        for (size_t i = 0; i < v; ++i) {
            T maxSimilarity = -numeric_limits<T>::max();
            int bestCluster = 0;
            for (size_t j = 0; j < k; ++j) {
                if (similarity[i][j] > maxSimilarity) {
                    maxSimilarity = similarity[i][j];
                    bestCluster = j;
                }
            }
            if (labels[i] != bestCluster) {
                labels[i] = bestCluster;
                changed = true;
            }
        }

        // 更新质心
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < d; j++) {
                centroids[i][j] = 0;
            }
        }
        vector<int> count(k, 0);
        for (size_t i = 0; i < v; ++i) {
            count[labels[i]]++;
            for (size_t j = 0; j < d; ++j) {
                centroids[labels[i]][j] += points[i][j];
            }
        }
        for (size_t j = 0; j < k; ++j) {
            if (count[j] > 0) {
                for (size_t l = 0; l < d; ++l) {
                    centroids[j][l] /= count[j];
                }
            }
        }
        Matrix<T>::normalizeRows(centroids);
        W_T = centroids.transpose();
        //if (!changed) break; //如果没有变化，提前终止
    }

    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "迭代次数： " << iter << endl;
    cout << endl;
    /*
    vector<queue<int>> queues(k);
    for (size_t i = 0; i < v; ++i) {
        queues[labels[i]].push(i);
    }
    for (size_t i = 0; i < k; i++) {
        cout << "簇" << i << "包含数据点：";
        while (!queues[i].empty()) {
            cout << queues[i].front() << " ";
            queues[i].pop();
        }
        cout << endl;
    }
    */
    cout << endl;
    cout << "k-means原始计算用时：" << t.count()*5 << endl;

    cout << endl;

    return { points, centroids };
}


//连乘原始计算
template <typename T>
Matrix<T> cMultiply_ori(Matrix<T>& X, Matrix<T>& A, Matrix<T>& B) {
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    chrono::high_resolution_clock::time_point start2;
    chrono::high_resolution_clock::time_point end2;
    chrono::milliseconds t;

    size_t n = X.getRows(), v = A.getRows(), d = A.getCols(), k = B.getRows();

    start = chrono::high_resolution_clock::now();
    //第一步：单位化
    Matrix<T>::normalizeRows(A);
    Matrix<T>::normalizeRows(B);
    Matrix<T> W_D = A;
    Matrix<T> W_T = B.transpose();

    //第二步：计算
    start2 = chrono::high_resolution_clock::now();
    Matrix<T> M = W_D * W_T;
    Matrix<T> Z = X * M;
    end2 = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
    cout << "连乘直接计算用时:" << " " << t.count() << endl;
    cout << endl;

    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
   
    cout << "连乘原始计算用时：" << t.count() << endl;
    cout << endl;

    return Z;
}

