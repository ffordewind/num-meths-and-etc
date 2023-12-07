#include <QuadMatrix.h>
#include <iostream>
using namespace std;
constexpr auto eps = 1e-3;
constexpr auto N = 3;

int main() {
    vector<vector<double>> exmpl = { {3, 1, 1}, {1, 5, 1}, {1, 1, 7} };
    vector<vector<double>> test0 = { {0, 2, 3}, {1, 2, 4}, {4, 5, 6} };
    vector<vector<double>> test1 = { {8, 1, 1}, {1, 10, 1}, {1, 1, 12} };
    vector<vector<double>> test2 = { {-8, 1, 1}, {1, -10, 1}, {1, 1, -12} };
    vector<vector<double>> test3 = { {-8, 9, 10}, {11, -10, 7}, {10, 11, -12} };
    vector<vector<double>> test4 = { {8, 7, 7}, {7, 10, 7}, {7, 7, 12} };



    //пример
    QuadMatrix A(exmpl);
    vector<double> b = { 5, 7, 9 };
    //тест 0
    QuadMatrix A0(test0);
    vector<double> b0 = { 13, 17, 32 };
    //тест 1
    QuadMatrix A1(test1);
    vector<double> b1 = { 10, 12, 14 };
    //тест 2
    QuadMatrix A2(test2);
    vector<double> b2 = { -10, -12, -14 };
    //тест 3
    QuadMatrix A3(test3);
    vector<double> b3 = { 10, 12, 14 };
    //тест 4
    QuadMatrix A4(test4);
    vector<double> b4 = { 10, 12, 14 };
    //тест 5
    QuadMatrix A5_0 = E(N);
    QuadMatrix A5_1(1, N);
    vector<double> b5;
    for (int i = 0; i < A5_0.dim; i++) {
        for (int j = 0; j < A5_0.dim; j++) {
            if (i < j) {
                A5_0.matrix[i][j] = -1;
                A5_1.matrix[i][j] = -1;
            }
        }
        b5.push_back(-1);
    }
    b5[N - 1] = 1;
    QuadMatrix A5 = A5_0 + A5_1 * (eps * 6);



    //"точные" решения систем из тестов
    vector<double> t0 = { 1, 2, 3 };
    vector<double> t1 = { 1, 1, 1 };
    vector<double> t2 = { 1.616379310345, 1.504310344828, 1.426724137931 };
    vector<double> t3 = { 1.4462540716612378, 1.1661237785016287, 1.1074918566775244 };
    vector<double> t4 = { -0.0227272727272727, 0.6590909090909091, 0.7954545454545455 };
    vector<double> t5_3_4 = { 0, 0, 0, 0.9940357852882704 };
    vector<double> t5_3_5 = { 0, 0, 0, 0, 0.9940357852882704 };
    vector<double> t5_3_6 = { 0, 0, 0, 0, 0, 0.9940357852882704 };
    vector<double> t5_6_4 = { 0, 0, 0, 0.9999940000359998 };
    vector<double> t5_6_5 = { 0, 0, 0, 0, 0.9999940000359998 };
    vector<double> t5_6_6 = { 0, 0, 0, 0, 0, 0.9999940000359998 };



    vector<double> x1 = LUP(A4, b4);
    vector<double> x2 = QR_HH(A4, b4);
    vector<double> x3 = SimpleIteration(A4, b4, t4);
    vector<double> x4 = Zeidel(A4, b4, t4);
    cout << setw(12 * N) << left << "Matrix A:" << setw(12) << right << "vector b:" << "\n\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(12) << left << A4.matrix[i][j];
        }
        cout << setw(12) << right << b4[i] << "\n\n\n";
    }
    cout << endl << endl << setw(16) << left << "LU(P):" << setw(16) << left << "Householder:" << setw(16) << left << "Iteration:" << setw(16) << left << "Zeidel:" << "\n\n";
    for (int i = 0; i < N; i++) {
        cout << setw(16) << left << x1[i] << setw(16) << left << x2[i] << setw(16) << left << x3[i] << setw(16) << left << x4[i] << endl;
    }
}