#include "Matrix.h"
#include "SparseMatrix.h"
int N = 500;

template <typename T>
chrono::microseconds calculateTime(const function<void()>& function)
{
    auto start = chrono::high_resolution_clock::now();
    function();
    auto end = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::microseconds>(end - start);
}

int main() {

    /*
    cout << "置换行：" << endl;
    Matrix<double> B = P1 * A;
    B.print();
    cout << endl;

    cout << "置换列：" << endl;
    Matrix<double> C = A * P2;
    C.print();
    cout << endl;

    Matrix<double> G = P1 * H1 * A * H2 * P2;
    cout << endl;
    */
    /*
    
    Matrix<double> A = Matrix<double>::random(N, N, 0, 10);
    Matrix<double> B = Matrix<double>::random(N, N, 0, 10);
    Matrix<double> I = Matrix<double>::RPMatrix(N);
    */

    /*
    auto start1 = chrono::high_resolution_clock::now();
    Matrix<double> C1 = A * B;
    auto end1 = chrono::high_resolution_clock::now();
    auto t1 = chrono::duration_cast<chrono::microseconds>(end1 - start1);
    cout << "直接计算用时：" << t1.count() << "e-3ms" << endl;
    C1.print();
    cout << endl;

    auto start2 = chrono::high_resolution_clock::now();
    Matrix<double> C2 = Matrix<double>::Coppersmith(A, B);
    auto end2 = chrono::high_resolution_clock::now();
    auto t2 = chrono::duration_cast<chrono::microseconds>(end2 - start2);
    cout << "Coppersmith方法用时：" << t2.count() << "e-3ms" << endl;
    //C2.print();
    cout << endl;

    //加密
    Matrix<double> Ap = P1 * H1 * A * H2 * P2;
    Matrix<double> Bp = P2T * H2 * B * H3 * P3;

    //计算
    Matrix<double> C3p = Matrix<double>::Coppersmith(Ap, Bp);

    //解密
    Matrix<double> C3 = H1 * P1T * C3p * P3T * H3;
    cout << "解密结果为：" << endl;
    C3.print();
    cout << endl;

    H2.print();
    cout << endl;
    Matrix<double> t = H2 * H2;
    t.print();
    */
    
    /*
    cout << "直接计算用时:" << " " << calculateTime<Matrix<double>>([&A, &H1, &P1]() { P1 * H1 * A ; }).count()/1000<<"ms" << endl;

    cout << "索引改进计算用时（左乘）:" << " " << calculateTime<Matrix<double>>([&A, &H1, &p1]() { Matrix<double>::ind_mul_left(H1, A, p1); }).count() / 1000 << "ms" << endl;

    SparseMatrix<double> A_s = A.toSparse();
    SparseMatrix<double> H1_s = H1.toSparse();
    SparseMatrix<double> P1_s = P1.toSparse();
    cout << "哈希表改进计算用时:" << " " << calculateTime<Matrix<double>>([&A_s, &H1_s, &P1_s]() { (P1_s * H1_s) * A_s; }).count()/1000<<"ms" << endl;
    */
    
    /*
    cout << "客户端承担更多加法：" << endl;
    Matrix<double> R = Matrix<double>::Coppersmith_1(A, B);

    cout << endl;

    cout << "服务器承担更多加法：" << endl;
    Matrix<double> R2 = Matrix<double>::Coppersmith_2(A, B);
    */
    /*
    Matrix<double> data({
        {1.0, 2.0},
        {-5.0, 8.0},
        {-1.5, 1.8},
        {4.5, 7.5},
        {1.8, -2.5},
        {4.5, -6.8},
        {-1.8, -1.5},
        {-14.8, -11.0}
        });
        */


    
    int n = 800;
    int v = 1200;
    int d = 300;
    int k = 1000;

    for (int i = 0; i < 8; i++) {

        cout << "n= " << n << " :" << endl;

        Matrix<double> X = Matrix<double>::random(n, v, 0, 10);
        Matrix<double> data = Matrix<double>::random(v, d, 0, 10);

        pair<Matrix<double>, Matrix<double>> result = kMeans(data, k, 5);

        Matrix<double> W_D = result.first;
        Matrix<double> W_T = result.second;

        Matrix<double> Z = cMultiply(X, W_D, W_T);

        pair<Matrix<double>, Matrix<double>> result_ori = kMeans_ori(data, k, 5);

        Matrix<double> W_D_ori = result_ori.first;
        Matrix<double> W_T_ori = result_ori.second;

        Matrix<double> Z_ori = cMultiply_ori(X, W_D, W_T);




        cout << "------------------------------------------------------------------------" << endl;

        n += 800;
        v += 1200;
    }
    


    /*
    W_D.print();
    cout << endl;
    W_D_ori.print();
    cout << endl;
    W_T.print();
    cout << endl;
    W_T_ori.print();
    cout << endl;
    Z.print();
    cout << endl;
    Z_ori.print();
    */
    

    /*
    Matrix<double> test = Matrix<double>::random(6, 6, 0, 10);
    Matrix<double> Q = Matrix<double>::householder(6, 1, 10);

    auto p = Matrix<double>::RPermutation(6);
    for (int i = 0; i < p.size(); i++) {
        cout << p[i] << " ";
    }
    cout << endl << endl;

    auto R = Matrix<double>::rightPMatrix(p);
    R.print();
    cout << endl;

    auto T1 = test * Q * R;
    T1.print();
    cout << endl;

    auto temp= Matrix<double>::ind_mul_r(test, Q);
    auto T2 =  temp * R;
    T2.print();
    cout << endl;

    auto T3= Matrix<double>::ind_mul_right(test, Q, p);
    T3.print();
    cout << endl;
    */



    return 0;
}
