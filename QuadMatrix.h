#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
constexpr auto e = 1e-2;

class QuadMatrix {
public:
	std::vector<std::vector<double>> matrix;
	int dim;

	QuadMatrix(std::vector<std::vector<double>> elems) { //матрица заданной размерности с заданными элементами
		this->dim = elems.size();
		for (int i = 0; i < dim; i++) {
			std::vector<double> temp;
			for (int j = 0; j < dim; j++) {
				temp.push_back(elems[i][j]);
			}
			this->matrix.push_back(temp);
		}
	};
	QuadMatrix(double e, int d) { //матрица заданной размерности с одинаковыми заданными элементами
		this->dim = d;
		for (int i = 0; i < d; i++) {
			std::vector<double> temp;
			for (int j = 0; j < d; j++) {
				temp.push_back(e);
			}
			this->matrix.push_back(temp);
		}
	};

	bool isSymmetric() { //проверка матрицы на симметричность
		bool check = 1;
		for (int i = 0; i < this->dim; i++) {
			for (int j = 0; j < this->dim; j++) {
				if (this->matrix[i][j] != this->matrix[j][i]) {
					check = 0;
				}
			}
		}
		return check;
	};
	double norm1() { //матричная норма 1
		double maximum = INT_MIN;
		for (int i = 0; i < this->dim; i++) {
			double sum = 0;
			for (int j = 0; j < this->dim; j++) {
				sum += abs(this->matrix[j][i]);
			}
			maximum = std::max(maximum, sum);
		}
		return maximum;
	};
	double normInf() { //матричная норма бесконечность
		double maximum = INT_MIN;
		for (int i = 0; i < this->dim; i++) {
			double sum = 0;
			for (int j = 0; j < this->dim; j++) {
				sum += abs(this->matrix[i][j]);
			}
			maximum = std::max(maximum, sum);
		}
		return maximum;
	};
	std::vector<double> getVector(int j, int k) { //взятие вектор-столбца из матрицы
		std::vector<double> v;
		for (int i = k; i < this->dim; i++) {
			v.push_back(this->matrix[i][j]);
		}
		return v;
	};
	QuadMatrix minor(int k) {
		QuadMatrix M(0, this->dim - k);
		for (int i = k; i < this->dim; i++) {
			for (int j = k; j < this->dim; j++) {
				M.matrix[i - k][j - k] = this->matrix[i][j];
			}
		}
		return M;
	};
	double det() { //определитель матрицы
		if (this->dim == 1)
			return this->matrix[0][0];
		else if (this->dim == 2)
			return this->matrix[0][0] * this->matrix[1][1] - this->matrix[0][1] * this->matrix[1][0];
		else {
			int d = 0;
			for (int k = 0; k < this->dim; k++) {
				QuadMatrix B(0, this->dim - 1);
				for (int i = 1; i < this->dim; i++) {
					for (int j = 0; j < this->dim; j++) {
						if (j == k) { continue; }
						else if (j < k) {
							B.matrix[i - 1][j] = this->matrix[i][j];
						}
						else {
							B.matrix[i - 1][j - 1] = this->matrix[i][j];
						}
					}
				}
				d += pow(-1, k + 2) * this->matrix[0][k] * B.det();
			}
			return d;
		}
	};
	bool diagDominant() { //проверка на диагональное преобладание
		for (int i = 0; i < this->dim; i++) {
			double sum = 0;
			for (int j = 0; j < this->dim; j++) {
				sum += this->matrix[i][j];
			}
			sum -= this->matrix[i][i];
			if (abs(sum) > abs(this->matrix[i][i])) { return 0; }
		}
		return 1;
	};

	QuadMatrix operator +(QuadMatrix B) { //сложение матриц
		for (int i = 0; i < this->dim; i++) {
			for (int j = 0; j < this->dim; j++) {
				B.matrix[i][j] += this->matrix[i][j];
			}
		}
		return B;
	};
	QuadMatrix operator -(QuadMatrix B) { //вычитание матриц
		for (int i = 0; i < this->dim; i++) {
			for (int j = 0; j < this->dim; j++) {
				B.matrix[i][j] -= this->matrix[i][j];
				if(B.matrix[i][j] != 0) { B.matrix[i][j] *= -1; }
			}
		}
		return B;
	};
	QuadMatrix operator *(QuadMatrix B) { //умножение матриц
		std::vector<std::vector<double>> prod;
		for (int i = 0; i < this->dim; i++) {
			std::vector<double> temp;
			for (int j = 0; j < this->dim; j++) {
				double sum = 0;
				for (int k = 0; k < this->dim; k++) {
					sum += this->matrix[i][k] * B.matrix[k][j];
				}
				temp.push_back(sum);
			}
			prod.push_back(temp);
		}
		QuadMatrix Prod(prod);
		return Prod;
	};
	std::vector<double> operator *(std::vector<double> X) { //умножение матрицы на вектор
		std::vector<double> prod;
		for (int i = 0; i < this->dim; i++) {
			double sum = 0;
			for (int j = 0; j < this->dim; j++) {
				sum += this->matrix[i][j] * X[j];
			}
			prod.push_back(sum);
		}
		return prod;
	};
	QuadMatrix operator *(double a) { //умножение матрицы на скаляр
		QuadMatrix Prod(this->matrix);
		for (int i = 0; i < this->dim; i++) {
			for (int j = 0; j < this->dim; j++) {
				Prod.matrix[i][j] *= a;
			}
		}
		return Prod;
	};
};

QuadMatrix E(int d) { //единичная матрица
	QuadMatrix I(0, d);
	for (int i = 0; i < I.dim; i++) {
		I.matrix[i][i] = 1;
	}
	return I;
}

QuadMatrix swapRows(QuadMatrix A, int i, int j) { //перестановка строк
	QuadMatrix As = A;
	As.matrix[i] = A.matrix[j];
	As.matrix[j] = A.matrix[i];
	return As;
}

QuadMatrix transpose(QuadMatrix A) { //транспонирование
	QuadMatrix At = A;
	for (int i = 0; i < A.dim; i++) {
		for (int j = 0; j < A.dim; j++) {
			At.matrix[i][j] = A.matrix[j][i];
			At.matrix[j][i] = A.matrix[i][j];
		}
	}
	return At;
}

double inProd(std::vector<double> A, std::vector<double> B) { //скалярное произведение
	double sum = 0;
	for (int i = 0; i < A.size(); i++) {
		sum += A[i] * B[i];
	}
	return sum;
}

QuadMatrix outProd(std::vector<double> A, std::vector<double> B) { //внешнее произведение векторов
	QuadMatrix Prod(0, A.size());
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A.size(); j++) {
			Prod.matrix[i][j] = A[i] * B[j];
		}
	}
	return Prod;
}

double root(double x) { //квадратный корень
	double s, s1, i = INT_MAX;
	s = (1 + x) / 2;
	while (i >= 1e-8) {
		s1 = (s + x / s) / 2;
		i = abs(s - s1);
		s = s1;
	}
	return s;
}

double normEuc_v(std::vector<double> A) { //евклидова векторная норма через квадратный корень выше
	double sum = 0;
	for (int i = 0; i < A.size(); i++) {
		sum += A[i] * A[i];
	}
	return root(sum);
}

double norm1_v(std::vector<double> A) { //векторная норма 1
	double sum = 0;
	for (int i = 0; i < A.size(); i++) {
		sum += abs(A[i]);
	}
	return sum;
}

double normInf_v(std::vector<double> A) { //векторная норма бесконечность
	double maximum = INT_MIN;
	for (int i = 0; i < A.size(); i++) {
		maximum = std::max(abs(A[i]), maximum);
	}
	return maximum;
}

std::vector<double> multiply(std::vector<double> a, double b) { //умножение вектора на скаляр
	std::vector<double> p;
	for (int i = 0; i < a.size(); i++) {
		p.push_back(a[i] * b);
	}
	return p;
}

std::vector<double> divide(std::vector<double> a, double b) { //деление вектора на скаляр
	std::vector<double> d;
	for (int i = 0; i < a.size(); i++) {
		d.push_back(a[i] / b);
	}
	return d;
}

std::vector<double> add_v(std::vector<double> a, std::vector<double> b) { //сумма векторов
	std::vector<double> c;
	for (int i = 0; i < a.size(); i++) {
		c.push_back(a[i] + b[i]);
	}
	return c;
}


std::vector<double> subtract_v(std::vector<double> a, std::vector<double> b) { //разность векторов
	std::vector<double> c;
	for (int i = 0; i < a.size(); i++) {
		c.push_back(a[i] - b[i]);
	}
	return c;
}










std::vector<double> LUP(QuadMatrix A, std::vector<double> b) { //LU(P) разложение матрицы и решение полученной системы
	QuadMatrix M(0, A.dim), P(0, A.dim);
	M = A;
	for (int i = 0; i < P.dim; i++) {
		P.matrix[i][i] = 1;
	}
	for (int i = 0; i < A.dim; i++) {
		double leadEl = 0;
		int leadRow = -1;
		for (int j = i; j < A.dim; j++) {
			if (fabs(M.matrix[j][i]) > leadEl) {
				leadEl = fabs(M.matrix[j][i]);
				leadRow = j;
			}
		}
		if (leadEl != 0) {
			P = swapRows(P, i, leadRow);
			M = swapRows(M, i, leadRow);
			for (int j = i + 1; j < A.dim; j++) {
				M.matrix[j][i] /= M.matrix[i][i];
				for (int k = i + 1; k < A.dim; k++) {
					M.matrix[j][k] -= M.matrix[j][i] * M.matrix[i][k];
				}
			}
		}
	}
	QuadMatrix U(0, A.dim), L(0, A.dim);
	for (int i = 0; i < M.dim; i++) {
		for (int j = 0; j < M.dim; j++) {
			if (i <= j) {
				U.matrix[i][j] = M.matrix[i][j];
			}
			else L.matrix[i][j] = M.matrix[i][j];
		}
		L.matrix[i][i] = 1;
	}
	std::vector<double> y(A.dim), x(A.dim);
	std::vector<double> c = P * b;
	for (int i = 0; i < L.dim; i++) {
		double sum = 0;
		for (int j = 0; j < i; j++) {
			sum += L.matrix[i][j] * y[j];
		}
		y[i] = (c[i] - sum) / L.matrix[i][i];
	}
	for (int i = U.dim - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < U.dim; j++) {
			sum += U.matrix[i][j] * x[j];
		}
		x[i] = (y[i] - sum) / U.matrix[i][i];
	}
	return x;
}



std::vector<double> QR_HH(QuadMatrix A, std::vector<double> b) { //QR разложение Хаусхолдера матрицы и решение полученной системы
	QuadMatrix Q = E(A.dim), Qtemp = E(A.dim), R = A, Rtemp = A;
	for (int i = 0; i < A.dim - 1; i++) {
		Qtemp = E(A.dim);
		std::vector<double> z(A.dim - i);
		z[0] = 1;
		std::vector<double> y1 = R.getVector(i, i);
		double a = normEuc_v(y1);
		double p = normEuc_v(subtract_v(y1, multiply(z, a)));
		std::vector<double> w = divide(subtract_v(y1, multiply(z, a)), p);
		QuadMatrix Q1 = E(A.dim - i) - outProd(w, w) * 2;
		Rtemp = Q1 * R.minor(i);
		for (int k = 0; k < R.dim - i; k++) {
			for (int l = 0; l < R.dim - i; l++) {
				R.matrix[k + i][l + i] = Rtemp.matrix[k][l];
				Qtemp.matrix[k + i][l + i] = Q1.matrix[k][l];
			}
		}
		Q = Q * Qtemp;
	}
	std::vector<double> y = transpose(Q) * b;
	std::vector<double> x(R.dim);
	for (int i = R.dim - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < R.dim; j++) {
			sum += R.matrix[i][j] * x[j];
		}
		x[i] = (y[i] - sum) / R.matrix[i][i];
	}
	return x;
}



double select_norm_v(int i, std::vector<double> v) { //выбор нормы вектора
	if (i == 0) { return norm1_v(v); }
	else return normInf_v(v);
}

std::vector<double> SimpleIteration(QuadMatrix A, std::vector<double> b, std::vector<double> truth) { //решение системы методом простой итерации
	int n;
	for (int i = 0; i < A.dim; i++) {
		if (A.minor(i).det() <= 0) {
			b = transpose(A) * b;
			A = transpose(A) * A;
		}
	}
	double norm = 0;
	double mu = 1 / A.normInf();
	n = 1;
	QuadMatrix B = E(A.dim) - A * mu;
	if (B.normInf() < 1) { norm = B.normInf(); }
	else {
		mu = 1 / A.norm1();
		n = 0;
		B = E(A.dim) - A * mu;
		if (B.norm1() < 1) { norm = B.norm1(); }
		else norm = 0.99;
	}
	std::vector<double> c = multiply(b, mu);
	std::vector<double> xold = c;
	std::cout << std::setw(24) << std::left << "Iteration iters" << std::setw(24) << std::left << "Iteration vectors:" << std::setw(24) << std::left << "Iteration deltas:" << "\n\n";
	for (int k = 0;;k++) {
		std::vector<double> x = add_v(B * xold, c);
		double crit1 = norm / (1 - norm) * select_norm_v(n, subtract_v(x, xold));
		double crit2 = 0;
		if (norm == 0.99) {
			crit2 = select_norm_v(n, subtract_v(A * x, b));
		}
		for (int i = 0; i < x.size(); i++) {
			if (i == 0) {
				std::cout << std::setw(24) << std::left << k << std::setw(24) << std::left << x[i] << std::setw(24) << std::left << abs(x[i] - truth[i]) << std::endl;
			}
			else std::cout << std::setw(24) << std::left << "" << std::setw(24) << std::left << x[i] << std::setw(24) << std::left << abs(x[i] - truth[i]) << std::endl;
		}
		std::cout << "\n\n\n";
		if (crit1 <= e && crit2 < e) {
			return x;
		}
		xold = x;
	}
}



std::vector<double> Zeidel(QuadMatrix A, std::vector<double> b, std::vector<double> truth) { //решение системы с помощью итераций Зейделя
	if (!A.diagDominant()) {
		b = transpose(A) * b;
		A = transpose(A) * A;
	}
	std::vector<double> d(A.dim);
	QuadMatrix C(0, A.dim);
	for (int i = 0; i < C.dim; i++) {
		d[i] = b[i] / A.matrix[i][i];
		for (int j = 0; j < C.dim; j++) {
			if (j != i) {
				C.matrix[i][j] = -A.matrix[i][j] / A.matrix[i][i];
			}
		}
	}
	std::vector<double> xold = d;
	std::cout << std::setw(24) << std::left << "Zeidel iters" << std::setw(24) << std::left << "Zeidel vectors:" << std::setw(24) << std::left << "Zeidel deltas:" << "\n\n";
	for (int k = 0;;k++) {
		std::vector<double> x(A.dim);
		for (int i = 0; i < A.dim; i++) {
			double sum = 0;
			for (int j = 0; j < A.dim; j++) {
				if (j >= i) {
					sum += C.matrix[i][j] * xold[j];
				}
				else sum += C.matrix[i][j] * x[j];
			}
			sum += d[i];
			x[i] = sum;
		}
		for (int i = 0; i < x.size(); i++) {
			if (i == 0) {
				std::cout << std::setw(24) << std::left << k << std::setw(24) << std::left << x[i] << std::setw(24) << std::left << abs(x[i] - truth[i]) << std::endl;
			}
			else std::cout << std::setw(24) << std::left << "" << std::setw(24) << std::left << x[i] << std::setw(24) << std::left << abs(x[i] - truth[i]) << std::endl;
		}
		std::cout << "\n\n\n";
		if (norm1_v(subtract_v(A * x, b)) <= e && normInf_v(subtract_v(A * x, b)) <= e) {
			return x;
		}
		xold = x;
	}
}