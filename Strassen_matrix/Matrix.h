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
class SparseMatrix; //ǰ������

template <class T>
class Matrix {
private:
    vector<vector<T>> data;
    size_t rows=0;
    size_t cols=0;

public:
    //��ʼ��
    //Ĭ�Ϲ��캯��
    Matrix() = default;
    //��ʽ���캯��
    Matrix(size_t r, size_t c) : rows(r), cols(c), data(r, vector<T>(c)) {}

    //������ʼ��
    Matrix(const vector<T>& rowVector) : rows(1), cols(rowVector.size()), data(1, rowVector) {}

    SparseMatrix<T> toSparse() const;

    //�б��ʼ����
    Matrix(initializer_list<initializer_list<T>> init) {
        rows = init.size();
        cols = init.begin()->size();
        data.resize(rows);
        auto data_ptr = data.begin();
        for (auto row = init.begin(); row != init.end(); ++row, ++data_ptr) {
            data_ptr->assign(row->begin(), row->end());
        }
    }

    //����˻�
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            cout << cols << "!=" << other.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < cols; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;//���ȿ���
            }
        }
        return result;
    }

    //if�Ż��˷�
    static Matrix<T> if_mul(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; ++i) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < a.cols; ++k) {
                    if (a.data[i][k] == 0 || b.data[k][j] == 0) result.data[i][j] = 0;
                    else result.data[i][j] += a.data[i][k] * b.data[k][j];
                }
                // ���ȿ���
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;
            }
        }
        return result;
    }

    //�����Ż��˷�(���)
    static Matrix<T> ind_mul_l(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i+=2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //�Կ���ÿ��
                    for (size_t l = 0; l < 2; ++l) {  //�Կ���ÿ��
                        result.data[i + k][j] += a.data[i + k][i + l] * b.data[i + l][j];
                    }
                    //if (abs(result.data[i + k][j]) < 1e-11) result.data[i + k][j] = 0;  //���ȿ���
                }
            }
        }
        return result;
    }

    //�����Ż��˷�(�ҳ�)
    static Matrix<T> ind_mul_r(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < result.cols; j += 2) {
                for (size_t k = 0; k < a.cols; ++k) {
                    for (size_t l = 0; l < 2; ++l) {
                        result.data[i][j + l] += a.data[i][k] * b.data[k][j + l];
                    }
                }
                //if (abs(result.data[i][j]) < 1e-11) result.data[i][j] = 0;  //���ȿ���
            }
        }
        return result;
    }

    //�����Ż��˷�(��ˣ��Ӵ��û����ȳ˺��û�)
    static Matrix<T> ind_mul_left(const Matrix& a, const Matrix& b, vector<int> p) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i += 2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //�Կ���ÿ��
                    for (size_t l = 0; l < 2; ++l) {  //�Կ���ÿ��
                        result.data[p[i + k]][j] += a.data[i + k][i + l] * b.data[i + l][j];
                    }
                }
                //if (abs(result.data[p[i]][j]) < 1e-11) result.data[p[i]][j] = 0;  //���ȿ���
            }
        }
        return result;
    }

    //�����Ż��˷�(��ˣ��Ӵ��û������û����)
    static Matrix<T> ind_mul_left_recover(const Matrix& a, const Matrix& b, vector<int> p) {
        if (a.cols != b.rows) {
            cout << a.cols << "!=" << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < result.rows; i += 2) {
            for (size_t j = 0; j < result.cols; ++j) {
                for (size_t k = 0; k < 2; ++k) {  //�Կ���ÿ��
                    for (size_t l = 0; l < 2; ++l) {  //�Կ���ÿ��
                        result.data[i + k][j] += a.data[i + k][i + l] * b.data[p[i + l]][j];
                    }
                }
                //if (abs(result.data[p[i]][j]) < 1e-11) result.data[p[i]][j] = 0;  //���ȿ���
            }
        }
        return result;
    }

    //�����Ż��˷�(�ҳˣ��Ӵ��û����ȳ˺��û�)
    static Matrix<T> ind_mul_right(const Matrix& a, const Matrix& b, const vector<int>& p) {
        // ������˷���ά���Ƿ�ƥ��
        if (a.cols != b.rows) {
            cout << a.cols << " != " << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols); 
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < b.cols; j += 2) {
                for (size_t k = 0; k < 2 ; ++k) { // �Կ���ÿ��
                    for (size_t l = 0; l < 2; ++l) { // �Կ���ÿ��
                        result.data[i][p[j + k]] += a.data[i][j + l] * b.data[j + l][j + k];
                    }
                    // ���ȿ���
                    // if (abs(result.data[i][p[j + k]]) < 1e-11) result.data[i][p[j + k]] = 0;
                }
            }
        }

        return result;
    }

    //�����Ż��˷�(�ҳˣ��Ӵ��û������û����)
    static Matrix<T> ind_mul_right_recover(const Matrix& a, const Matrix& b, const vector<int>& p) {
        // ������˷���ά���Ƿ�ƥ��
        if (a.cols != b.rows) {
            cout << a.cols << " != " << b.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }
        Matrix result(a.rows, b.cols);
        for (size_t i = 0; i < a.rows; ++i) {
            for (size_t j = 0; j < b.cols; j += 2) {
                for (size_t k = 0; k < 2; ++k) { // �Կ���ÿ��
                    for (size_t l = 0; l < 2; ++l) { // �Կ���ÿ��
                        result.data[i][j + k] += a.data[i][p[j + l]] * b.data[j + l][j + k];
                    }
                    // ���ȿ���
                    // if (abs(result.data[i][p[j + k]]) < 1e-11) result.data[i][p[j + k]] = 0;
                }
            }
        }

        return result;
    }
  

    Matrix<T> mulByr(const Matrix<T>& other) const {
        if (cols != other.getRows() || other.getCols() != 1) {
            cout << cols << "!=" << other.rows << endl;
            throw invalid_argument("�˷������в�ƥ��");
        }

        Matrix<T> result(rows, 1);
        for (size_t j = 0; j < cols; ++j) {
            if (other.data[j][0] == 1) {
                for (size_t i = 0; i < rows; ++i) {
                    result.data[i][0] += data[i][j];  //�ۼ�A�ĵ�j�е��������
                }
            }
        }
        return result;
    }

    //��������
    Matrix operator*(const T& scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    //����ӷ�
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("�ӷ����в�ƥ��");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    //�������
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("�������в�ƥ��");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    //����ת��
    Matrix<T> transpose() const {
        Matrix<T> result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    //�������
    Matrix<T> inverse() const {
        if (rows != cols) {
            throw invalid_argument("��������Ƿ���");
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

        for (size_t i = 0; i < rows; ++i) { //��λ����
            result[i][i] = static_cast<T>(1);
        }

        for (size_t i = 0; i < rows; ++i) {
            // Ѱ����Ԫ
            if (copy[i][i] == static_cast<T>(0)) {
                bool found = false;
                for (size_t j = i + 1; j < rows && !found; ++j) {
                    if (copy[j][i] != static_cast<T>(0)) {
                        // ������
                        swap(copy.data[i], copy.data[j]);
                        swap(result.data[i], result.data[j]);
                        found = true;
                    }
                }
                if (!found) {
                    throw runtime_error("���󲻿���");
                }
            }

            // ʹ��Ԫ��Ϊ1
            T div = copy[i][i];
            for (size_t j = 0; j < cols; ++j) {
                copy[i][j] /= div;
                result[i][j] /= div;
            }

            // ��ȥ�����еĸ���Ԫ��
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

    //�������
    vector<T>& operator[](size_t index) {
        if (index >= rows) {
            throw std::out_of_range("����Խ��");
        }
        return data[index];
    }

    //�����ֵ�Ƚ�
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) return false;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (abs(data[i][j] - other.data[i][j]) > 1e-3) return false;
            }
        }
        return true;
    }

    //����ڶ���ʽ��λ��(||A||)
    static void normalizeRows(Matrix<T>& matrix) {
        for (size_t i = 0; i < matrix.getRows(); ++i) {
            T norm = 0;
            for (size_t j = 0; j < matrix.getCols(); ++j) {
                norm += matrix[i][j] * matrix[i][j];
            }
            norm = sqrt(norm);
            if (norm > 0) { // ��ֹ������
                for (size_t j = 0; j < matrix.getCols(); ++j) {
                    matrix[i][j] /= norm;
                }
            }
        }
    }

    //��ӡ
    void print() const {
        for (const auto& row : data) {
            for (const auto& elem : row) {
                cout << elem << " ";
            }
            cout << endl;
        }
    }

    //��ȡ����
    size_t getRows() const { return rows; }

    //��ȡ����
    size_t getCols() const { return cols; }

    //��ȡ������Ϊ����
    Matrix<T> getRow(size_t index) const {
        if (index >= rows) {
            throw std::out_of_range("����Խ��");
        }
        return Matrix<T>{data[index]};
    }

    //ŷ�Ͼ���
    static T distance(Matrix a, Matrix b) {
        if (a.getCols() != b.getCols() || a.getRows() != 1 || b.getRows() != 1) {
            throw std::invalid_argument("��������ƥ��");
        }
        T dist = 0;
        for (size_t i = 0; i < a.getCols(); i++) {
            T diff = a[0][i] - b[0][i];
            dist += diff * diff;
        }
        return sqrt(dist);
    }

    //�������ƶ�
    static Matrix<T> cos(Matrix<T> A, Matrix<T> B) {
        Matrix Anorm = A, Bnorm = B;
        //��λ��
        normalizeRows(Anorm);
        normalizeRows(Bnorm);

        return Anorm * Bnorm.transpose();
    }

    //�����������
    static Matrix random(size_t r, size_t c, T min = numeric_limits<T>::min(), T max = numeric_limits<T>::max()) {
        random_device rd;
        mt19937 gen(rd()); //α�����������
        uniform_real_distribution<> dis(min, max);

        Matrix randomMatrix(r, c);
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < c; ++j) {
                //�����������������Ҫת��
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

    //Householder�任
    static Matrix<T> householder(size_t n, T min, T max) {
        Matrix<T> result(n, n);
        if (n % 2 != 0)  throw runtime_error("������2�ı���");
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

    //��������û�����
    static Matrix<T> RPMatrix(size_t size) {
        Matrix<T> PMatrix(size, size);

        //��䵥λ����
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[i][i] = static_cast<T>(1);
        }

        //�����������
        random_device rd;
        mt19937 generator(rd());

        vector<int> indices(size);
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), generator);
        //�û�
        for (size_t i = 0; i < size; ++i) {
            if (i != indices[i]) {
                swap(PMatrix.data[i], PMatrix.data[indices[i]]);
            }
        }
        return PMatrix;
    }

    //����û�ӳ��
    static vector<int> RPermutation(size_t size) {
        vector<int> i(size);
        iota(i.begin(), i.end(), 0);
        random_device rd;
        mt19937 generator(rd());
        shuffle(i.begin(), i.end(), generator);

        return i;
    }

    //����û�������ˣ�������(P*A)
    static Matrix<T> leftPMatrix(const vector<int>& p) {
        size_t size = p.size();
        Matrix<T> PMatrix(size, size);
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[p[i]][i] = static_cast<T>(1);  //��1�����û�ָ����
        }
        return PMatrix;
    }

    //����û������ҳˣ�������(A*P)
    static Matrix<T> rightPMatrix(const vector<int>& p) {
        size_t size = p.size();
        Matrix<T> PMatrix(size, size);
        for (size_t i = 0; i < size; ++i) {
            PMatrix.data[i][p[i]] = static_cast<T>(1);  //��1�����û�ָ����
        }
        return PMatrix;
    }

    //���������֤����
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


    //strassen�㷨
    static Matrix<T> Strassen(Matrix<T> A, Matrix<T> B, int &ss, int &cc) {
        if (A.getCols() % 2 != 0 || B.getCols() % 2 != 0) throw runtime_error("������2�ı���");
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
        // << "�ͻ��˼ӷ�Ԥ���������ʱ:" << " " << t1.count() << endl;
        ss += t1.count();
        //��������ɵĲ���----------------------------------------
        auto start2 = chrono::high_resolution_clock::now();
        Matrix<T> M1 = A11 * R1;
        Matrix<T> M2 = R2 * B22;
        auto end2 = chrono::high_resolution_clock::now();
        auto t2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        //cout << "S1������ʱ:" << " " << t2.count() << endl;

        auto start3 = chrono::high_resolution_clock::now();
        Matrix<T> M3 = R3 * B11;
        Matrix<T> M4 = A22 * R4;
        auto end3 = chrono::high_resolution_clock::now();
        auto t3 = chrono::duration_cast<chrono::milliseconds>(end3 - start3);
        //cout << "S2������ʱ:" << " " << t3.count() << endl;

        auto start4 = chrono::high_resolution_clock::now();
        Matrix<T> M5 = R5 * R6;
        Matrix<T> M6 = R7 * R8;
        auto end4 = chrono::high_resolution_clock::now();
        auto t4 = chrono::duration_cast<chrono::milliseconds>(end4 - start4);
        //cout << "S3������ʱ:" << " " << t4.count() << endl;

        auto start5 = chrono::high_resolution_clock::now();
        Matrix<T> M7 = R9* R10;    
        auto end5 = chrono::high_resolution_clock::now();
        auto t5 = chrono::duration_cast<chrono::milliseconds>(end5 - start5);
        //cout << "S4������ʱ:" << " " << t5.count() << endl;

        auto tt = t2 + t3 + t4 + t5 ;
        //cout << "strassen�����Ʒ���������ʱ:" << " " << tt.count() << endl;
        ss += tt.count();

        //--------------------------------------------------------
        auto start9 = chrono::high_resolution_clock::now();
        Matrix<T> U11 = M5 + M4 - M2 + M6;
        Matrix<T> U12 = M1 + M2;
        Matrix<T> U21 = M3 + M4;
        Matrix<T> U22 = M5 + M1 - M3 - M7;
        auto end9 = chrono::high_resolution_clock::now();
        auto t9 = chrono::duration_cast<chrono::milliseconds>(end9 - start9);
        //cout << "�ͻ��˼ӷ����������ʱ:" << " " << t9.count() << endl;
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

    //Coppersmith-Winograd�㷨(һ�֣��ͻ��˳е�����ӷ�)
    static Matrix<T> Coppersmith(Matrix<T> A, Matrix<T> B, int &s, int &c) {
        if (A.getCols() % 2 != 0 || B.getCols() % 2 != 0) throw runtime_error("������2�ı���");
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
        // << "�ͻ��˼ӷ�Ԥ���������ʱ:" << " " << t1.count() << endl;
        c += t1.count();
        //��������ɵĲ���----------------------------------------
        auto start2 = chrono::high_resolution_clock::now();
        Matrix<T> M1 = A11 * B11;
        auto end2 = chrono::high_resolution_clock::now();
        auto t2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        //cout << "S1������ʱ:" << " " << t2.count() << endl;

        auto start3 = chrono::high_resolution_clock::now();
        Matrix<T> M2 = A12 * B21;
        auto end3 = chrono::high_resolution_clock::now();
        auto t3 = chrono::duration_cast<chrono::milliseconds>(end3 - start3);
        //cout << "S2������ʱ:" << " " << t3.count() << endl;

        auto start4 = chrono::high_resolution_clock::now();
        Matrix<T> S4 = A12 - S2;
        Matrix<T> M3 = S4 * B22;
        auto end4 = chrono::high_resolution_clock::now();
        auto t4 = chrono::duration_cast<chrono::milliseconds>(end4 - start4);
        //cout << "S3������ʱ:" << " " << t4.count() << endl;

        auto start5 = chrono::high_resolution_clock::now();
        Matrix<T> T4 = T2 - B21;
        Matrix<T> M4 = A22 * T4;
        auto end5 = chrono::high_resolution_clock::now();
        auto t5 = chrono::duration_cast<chrono::milliseconds>(end5 - start5);
        //cout << "S4������ʱ:" << " " << t5.count() << endl;

        auto start6 = chrono::high_resolution_clock::now();
        Matrix<T> M5 = S1 * T1;
        auto end6 = chrono::high_resolution_clock::now();
        auto t6 = chrono::duration_cast<chrono::milliseconds>(end6 - start6);
        //cout << "S5������ʱ:" << " " << t6.count() << endl;

        auto start7 = chrono::high_resolution_clock::now();
        Matrix<T> M6 = S2 * T2;
        auto end7 = chrono::high_resolution_clock::now();
        auto t7 = chrono::duration_cast<chrono::milliseconds>(end7 - start7);
        //cout << "S6������ʱ:" << " " << t7.count() << endl;

        auto start8 = chrono::high_resolution_clock::now();
        Matrix<T> S3 = A11 - A21;
        Matrix<T> T3 = B22 - B12;
        Matrix<T> M7 = S3 * T3;
        auto end8 = chrono::high_resolution_clock::now();
        auto t8 = chrono::duration_cast<chrono::milliseconds>(end8 - start8);
        //cout << "S7������ʱ:" << " " << t8.count() << endl;

        auto tt = t2 + t3 + t4 + t5 + t6 + t7 + t8;
        //cout << "�Ʒ���������ʱ:" << " " << tt.count() << endl;
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
        //cout << "�ͻ��˼ӷ����������ʱ:" << " " << t9.count() << endl;
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

//����ͨ����ת��Ϊϡ�����
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















//K-means���లȫ���
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
    Matrix<T> centroids(k, d); //���ľ���
    vector<int> labels(v, 0);
    int sp = 10;

    //��һ������ʼ�������ѡ������
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, v - 1);
    for (size_t i = 0; i < k; ++i) {
        size_t randIdx = dis(gen);
        for (size_t j = 0; j < d; ++j) {
            centroids[i][j] = points[randIdx][j];
        }
    }

    //�ڶ�������λ��
    Matrix<T>::normalizeRows(points);
    Matrix<T>::normalizeRows(centroids);
    Matrix<T> W_D = points;
    Matrix<T> W_T = centroids.transpose();

    int iter = 0;
    for (; iter < 1; ++iter) {
        bool changed = false;

        //����������Կ����
        Matrix<double> Q1 = Matrix<double>::householder(v, -pow(2, sp), pow(2, sp));
        Matrix<double> Q2 = Matrix<double>::householder(d, -pow(2, sp), pow(2, sp));
        Matrix<double> Q3 = Matrix<double>::householder(k, -pow(2, sp), pow(2, sp));
        auto p1 = Matrix<double>::RPermutation(v);
        auto p2 = Matrix<double>::RPermutation(d);
        auto p3 = Matrix<double>::RPermutation(k);

        //���Ĳ���ä��
        start = chrono::high_resolution_clock::now();
        Matrix<T> W_D_t = Matrix<double>::ind_mul_left(Q1, W_D, p1);
        Matrix<T> W_T_t = Matrix<double>::ind_mul_left(Q2, W_T, p2);
        Matrix<T> W_D_b = Matrix<double>::ind_mul_right(W_D_t, Q2, p2);
        Matrix<T> W_T_b = Matrix<double>::ind_mul_right(W_T_t, Q3, p3);
        Matrix<T> r[10];
        Matrix<double> V[10];
        for (int i = 0; i < 10; i++) {
            r[i] = Matrix<T>::R01r(k);
            V[i] = W_T_b.mulByr(r[i]);//������֤
        }
        end = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end - start);
        blitime += t;
        cout << "ä����ʱ:" << " " << t.count() << endl;

        //���岽������
        Matrix<T> C_b = Matrix<double>::Coppersmith(W_D_b, W_T_b, s, c);

        Matrix<T> t1 = Matrix<double>::Strassen(W_D_b, W_T_b, ss, cc);

        //����������֤  
        t = chrono::duration_cast<chrono::milliseconds>(start - start);
        Matrix<double> temp,temp2;
        for (int i = 0; i < 10; i++) {
            start = chrono::high_resolution_clock::now();
            temp = C_b.mulByr(r[i]);
            temp2 = W_D_b * V[i];
            end = chrono::high_resolution_clock::now();
            t += chrono::duration_cast<chrono::milliseconds>(end - start);
            if (!(temp == temp2)) {
                cout << "��" << i + 1 << "����֤δͨ��" << endl;
                break;
            }
            else cout << "��" << i + 1 << "����֤ͨ��" << endl;

        }
        vertime += t;
        cout << "��֤��ʱ:" << " " << t.count() << endl;
       
        //���߲����ָ�
        start = chrono::high_resolution_clock::now();
        Matrix<T> C_h = Matrix<double>::ind_mul_left_recover(Q1, C_b, p1);
        Matrix<T> C = Matrix<double>::ind_mul_right_recover(C_h, Q3, p3);
        end = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end - start);
        rectime += t;
        cout << "�ָ���ʱ:" << " " << t.count() << endl;
        cout << endl;

        //�ڰ˲�������
        Matrix<T> similarity = C;
        // ���·��������Ƶ�����
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

        // ��������
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
        //if (!changed) break; //���û�б仯����ǰ��ֹ


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

    cout << "���������� " << iter*5  << endl;
    cout << endl;
    /*
    vector<queue<int>> queues(k);
    for (size_t i = 0; i < v; ++i) {
        queues[labels[i]].push(i);
    }
    for (size_t i = 0; i < k; i++) {
        cout << "��" << i << "�������ݵ㣺";
        while (!queues[i].empty()) {
            cout << queues[i].front() << " ";
            queues[i].pop();
        }
        cout << endl;
    }
    */
    cout << endl;
    cout << "ä��k-means����ʱ��" << blitime.count()*5 << endl;
    cout << "��֤k-means����ʱ��" << vertime.count()*5 << endl;
    cout << "�ָ�k-means����ʱ��" << rectime.count()*5 << endl;
    cout << "�Ʒ�����k-means����ʱ��" << s * 5 << endl;
    cout << "�ͻ���k-meansԤ������ʱ��" << c * 5 << endl;
    cout << endl;
    cout << "strassen��k-means����ʱ��" << ss * 5 << endl;
    cout << endl;
    cout << "liu����k-means����ʱ��" << liutime * 5 << endl;
    cout << endl;

    return { points, centroids };
}


//���˰�ȫ���
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

    //��һ������λ��
    Matrix<T>::normalizeRows(A);
    Matrix<T>::normalizeRows(B);
    Matrix<T> W_D = A;
    Matrix<T> W_T = B.transpose();

    //�ڶ�������Կ����
    Matrix<double> Q1 = Matrix<double>::householder(n, -pow(2, sp), pow(2, sp));
    Matrix<double> Q2 = Matrix<double>::householder(v, -pow(2, sp), pow(2, sp));
    Matrix<double> Q3 = Matrix<double>::householder(d, -pow(2, sp), pow(2, sp));
    Matrix<double> Q4 = Matrix<double>::householder(k, -pow(2, sp), pow(2, sp));
    auto p1 = Matrix<double>::RPermutation(n);
    auto p2 = Matrix<double>::RPermutation(v);
    auto p3 = Matrix<double>::RPermutation(d);
    auto p4 = Matrix<double>::RPermutation(k);

    //��������ä��
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
        V[i] = W_T_b.mulByr(r[i]);//������֤
    }
    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    liutime += t.count();
    cout << "ä����ʱ:" << " " << t.count() << endl;

    //���Ĳ�������
    Matrix<T> M = Matrix<double>::Coppersmith(W_D_b, W_T_b, s ,c);
    Matrix<T> Z_b = Matrix<double>::Coppersmith(X_b, M, s, c);

    Matrix<T> t1 = Matrix<double>::Strassen(W_D_b, W_T_b, ss, cc);
    Matrix<T> t2 = Matrix<double>::Strassen(X_b, t1, ss, cc);

    //���岽����֤
    t = chrono::duration_cast<chrono::milliseconds>(start - start);
    Matrix<double> temp, temp2;
    for (int i = 0; i < 10; i++) {
        start = chrono::high_resolution_clock::now();
        temp = Z_b.mulByr(r[i]);
        temp2 = X_b * (W_D_b * V[i]);
        end = chrono::high_resolution_clock::now();
        t += chrono::duration_cast<chrono::milliseconds>(end - start);
        if (!(temp == temp2)) {
            cout << "��" << i + 1 << "����֤δͨ��" << endl;
            break;
        }
        else cout << "��" << i + 1 << "����֤ͨ��" << endl;

    }
    liutime += t.count();
    cout << "��֤��ʱ:" << " " << t.count() << endl;

    //���������ָ�
    start = chrono::high_resolution_clock::now();
    Matrix<T> Z_h = Matrix<double>::ind_mul_left_recover(Q1, Z_b, p1);
    Matrix<T> Z = Matrix<double>::ind_mul_right_recover(Z_h, Q4, p4);
    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    liutime += t.count();
    cout << "�ָ���ʱ:" << " " << t.count() << endl;
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


    cout << "�Ʒ�����������ʱ��" << s << endl;
    cout << "�ͻ���Ԥ������ʱ��" << c << endl;
    cout << endl;
    cout << "strassen������ʱ��" << ss << endl;
    cout << endl;
    cout << "liu������������ʱ��" << liutime << endl;
    cout << endl;

    return Z;
}






//k-meansԭʼ����
template <typename T>
pair<Matrix<T>, Matrix<T>> kMeans_ori(Matrix<T>& points, size_t k, int maxIterations) {
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    chrono::high_resolution_clock::time_point start2;
    chrono::high_resolution_clock::time_point end2;
    chrono::milliseconds t;

    size_t v = points.getRows(), d = points.getCols();
    Matrix<T> centroids(k, d); //���ľ���
    vector<int> labels(v, 0);

    start = chrono::high_resolution_clock::now();
    //��һ������ʼ�������ѡ������
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, v - 1);
    for (size_t i = 0; i < k; ++i) {
        size_t randIdx = dis(gen);
        for (size_t j = 0; j < d; ++j) {
            centroids[i][j] = points[randIdx][j];
        }
    }

    //�ڶ�������λ��
    Matrix<T>::normalizeRows(points);
    Matrix<T>::normalizeRows(centroids);
    Matrix<T> W_D = points;
    Matrix<T> W_T = centroids.transpose();

    int iter = 0;
    for (; iter < 1; ++iter) {
        bool changed = false;

        //������������
        start2 = chrono::high_resolution_clock::now();
        Matrix<T> C = W_D * W_T;
        end2 = chrono::high_resolution_clock::now();
        t = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
        cout << "k-meansֱ�Ӽ�����ʱ:" << " " << t.count() * 5 << endl;
        cout << endl;

        //���Ĳ�������
        Matrix<T> similarity = C;
        // ���·��������Ƶ�����
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

        // ��������
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
        //if (!changed) break; //���û�б仯����ǰ��ֹ
    }

    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "���������� " << iter << endl;
    cout << endl;
    /*
    vector<queue<int>> queues(k);
    for (size_t i = 0; i < v; ++i) {
        queues[labels[i]].push(i);
    }
    for (size_t i = 0; i < k; i++) {
        cout << "��" << i << "�������ݵ㣺";
        while (!queues[i].empty()) {
            cout << queues[i].front() << " ";
            queues[i].pop();
        }
        cout << endl;
    }
    */
    cout << endl;
    cout << "k-meansԭʼ������ʱ��" << t.count()*5 << endl;

    cout << endl;

    return { points, centroids };
}


//����ԭʼ����
template <typename T>
Matrix<T> cMultiply_ori(Matrix<T>& X, Matrix<T>& A, Matrix<T>& B) {
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    chrono::high_resolution_clock::time_point start2;
    chrono::high_resolution_clock::time_point end2;
    chrono::milliseconds t;

    size_t n = X.getRows(), v = A.getRows(), d = A.getCols(), k = B.getRows();

    start = chrono::high_resolution_clock::now();
    //��һ������λ��
    Matrix<T>::normalizeRows(A);
    Matrix<T>::normalizeRows(B);
    Matrix<T> W_D = A;
    Matrix<T> W_T = B.transpose();

    //�ڶ���������
    start2 = chrono::high_resolution_clock::now();
    Matrix<T> M = W_D * W_T;
    Matrix<T> Z = X * M;
    end2 = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
    cout << "����ֱ�Ӽ�����ʱ:" << " " << t.count() << endl;
    cout << endl;

    end = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::milliseconds>(end - start);
   
    cout << "����ԭʼ������ʱ��" << t.count() << endl;
    cout << endl;

    return Z;
}

