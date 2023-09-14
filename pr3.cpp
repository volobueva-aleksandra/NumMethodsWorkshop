#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )
#define sign(x) ( (x) == 0 ? 0 : ( (x) < 0 ? (-1) : (1) ) )

#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0
#define _USE_MATH_DEFINES

#define n 13
#define xend 1

static double Runge_Kutta(void(*f)(double, double, double*, double*), double y[][n], double x,
	double h, double a);

void f(double a, double x, double* y, double* ans)
{
	ans[0] = y[1];
	ans[1] = y[3] * (1 + a * y[0] * y[0]);
	ans[2] = -y[3] * y[3] * a * y[0];
	ans[3] = -y[2];

	ans[4] = y[5];
	ans[5] = y[7] * (1 + a * y[0] * y[0]) + y[3] * 2 * a * y[0] * y[4];
	ans[6] = -2 * y[7] * a * y[0] * y[3] - y[3] * y[3] * a * y[4];
	ans[7] = -y[6];

	ans[8] = y[3] * y[3] *(1 + a * y[0] * y[0]);

	ans[9] = y[10];
	ans[10] = y[12] * (1 + a * y[0] * y[0]) + y[3] * 2 * a * y[0] * y[9];
	ans[11] = -2 * y[12] * a * y[0] * y[3] - y[3] * y[3] * a * y[9];
	ans[12] = -y[11];
}

double Prince_Dormand(void(*f)(double, double, double*, double*), double y[][n],
	double x, double h, double xmax, double* h_next, double tolerance, double a) {

	double scale;
	double temp_y[2][n];
	double err = 0;
	double yy = 0;
	int i, j, k = 0, g=-2;
	int last_interval = 0;

	if (xmax < x || h <= 0.0) return -2;
	*h_next = h;
	for (i = 0; i < n; i++)
		y[1][i] = y[0][i];
	if (xmax == x) return 0;
	h = min(h, xmax - x);
	tolerance /= (xmax - x);
	for (i = 0; i < n; i++)
		temp_y[0][i] = y[0][i];

	while (x < xmax) {
		scale = 1.0;
		for (i = 0; i < ATTEMPTS; i++) {
			yy = 0;
			err = fabs(Runge_Kutta(f, temp_y, x, h, a));
			if (err == 0.0) { scale = MAX_SCALE_FACTOR; break; }
			for (j = 0; j < n; j++) yy += (temp_y[0][j] == 0.0) ? tolerance : fabs(temp_y[0][j]);
			scale = 0.8 * sqrt(sqrt(tolerance * yy / err));
			scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
			if (err < (tolerance * yy)) break;
			h *= scale;
			if (x + h > xmax) h = xmax - x;
			else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
		}
		if (i >= ATTEMPTS) { *h_next = h * scale; return -1; };
		for (j = 0; j < n; j++) temp_y[0][j] = temp_y[1][j];
		x += h;
		h *= scale;
		*h_next = h;
		if (last_interval) break;
		if (x + h > xmax) { last_interval = 1; h = xmax - x; }
		else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
		k++;
		if (k == 2 and fabs(y[0][5]) > 0)
			g = sign(temp_y[1][9] * temp_y[1][5] - temp_y[1][10] * temp_y[1][4]);
		if (k > 2 and fabs(y[0][5]) > 0 && g * sign(temp_y[1][9] * temp_y[1][5] - temp_y[1][10] * temp_y[1][4]) <= 0)
			cout << "focal point = " << x << endl;
	}
	for (j = 0; j < n; j++) y[1][j] = temp_y[1][j];
	return err;
}


static double Runge_Kutta(void(*f)(double, double, double*, double*), double y[][n], double x0,
	double h, double a) {
	static const double r_45 = 1.0 / 45.0;
	static const double r_8_9 = 8.0 / 9.0;
	static const double r_6561 = 1.0 / 6561.0;
	static const double r_167904 = 1.0 / 167904.0;
	static const double r_142464 = 1.0 / 142464.0;
	static const double r_21369600 = 1.0 / 21369600.0;

	double y_tmp[n];
	double err = 0;

	double k1[n], k2[n], k3[n], k4[n], k5[n], k6[n], k7[n];
	double h5 = 0.2 * h;

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i];
	(*f)(a, x0, y_tmp, k1);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + h5 * k1[i];
	(*f)(a, x0 + h5, y_tmp, k2);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + h * (0.075 * k1[i] + 0.225 * k2[i]);
	(*f)(a, x0 + 0.3 * h, y_tmp, k3);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + h * r_45 * (44.0 * k1[i] - 168.0 * k2[i] + 160 * k3[i]);
	(*f)(a, x0 + 0.8 * h, y_tmp, k4);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + r_6561 * h * (19372.0 * k1[i]
			- 76080.0 * k2[i] + 64448.0 * k3[i] - 1908.0 * k4[i]);
	(*f)(a, x0 + r_8_9 * h, y_tmp, k5);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + r_167904 * h * (477901.0 * k1[i] - 1806240.0 * k2[i]
			+ 1495424.0 * k3[i] + 46746.0 * k4[i] - 45927.0 * k5[i]);
	(*f)(a, x0 + h, y_tmp, k6);

	for (int i = 0; i < n; i++)
		y_tmp[i] = y[0][i] + r_142464 * h * (12985.0 * k1[i] + 64000.0 * k3[i]
			+ 92750.0 * k4[i] - 45927.0 * k5[i] + 18656.0 * k6[i]);
	(*f)(a, x0 + h, y_tmp, k7);


	for (int i = 0; i < n; i++)
	{
		y[1][i] = y[0][i] + r_21369600 * h * (1921409.0 * k1[i] + 9690880.0 * k3[i]
			+ 13122270.0 * k4[i] - 5802111.0 * k5[i] + 1902912.0 * k6[i] + 534240.0 * k7[i]);
		err += fabs(r_21369600 * (26341.0 * k1[i] - 90880.0 * k3[i] + 790230.0 * k4[i]
			- 1086939.0 * k5[i] + 895488.0 * k6[i] - 534240.0 * k7[i]));
	}
	return err;
}

void shot(double a, double eps, double y[][n], double h, double alpha, double beta)
{
	double h_next;
	double err;
	for (int i = 0; i < n; i++)
		y[0][i] = 0;
	y[0][1] = 1; y[0][2] = alpha; y[0][3] = beta;
	err = Prince_Dormand(&f, y, 0, h, xend, &h_next, eps, a);
	(void)err;
}

double shot1(double a, double eps, double y[][n], double h, double alpha, double beta)
{
	double h_next;
	double err;
	for (int i = 0; i < n; i++)
		y[0][i] = 0;
	y[0][2] = alpha; y[0][3] = beta;
	y[0][1] = 1; y[0][6] = 1; y[0][12] = 1;
	err = Prince_Dormand(&f, y, 0, h, xend, &h_next, eps, a);
	return err;
}

int main()
{
	cout << " \\begin{table}[H] \\begin{center} \\begin{tabular}{| c || c | c | c |}" << endl;
	cout << " \\hline $\\alpha$ & $\\beta_1$ & $\\beta_2$ & $B_0$ \\\\ \\hline" << endl;
	double a = 0.1;
		double eps = pow(10, -5);
		double y[2][n];
		double h = 0.0001;
		double S = 1;
		double alpha = 0.5, beta = -0.5;
		while (S > eps)
		{
			double a1, a2, a3, a4, b1, b2;
			shot(a, eps, y, h, alpha, beta);

			b1 = -y[1][1]; b2 = -y[1][2];
			S = b1 * b1 + b2 * b2;
			shot(a, eps, y, h, alpha + eps, beta);
			a1 = (y[1][1] + b1) / eps; a3 = (y[1][2] + b2) / eps;
			shot(a, eps, y, h, alpha, beta + eps);
			a2 = (y[1][1] + b1) / eps; a4 = (y[1][2] + b2) / eps;
			if (fabs(a1 * a4 - a3 * a2) < eps)
			{
				cout << "Error";
				break;
			}
			alpha += (b1 * a4 - b2 * a2) / (a1 * a4 - a3 * a2);
			beta += (b2 * a1 - b1 * a3) / (a1 * a4 - a3 * a2);
		}
		double err = shot1(a, eps, y, h, alpha, beta);
		cout.precision(1);
		cout << fixed << a << "& ";
		cout.precision(10);
		cout << alpha << "& " << beta << "& " << y[1][8] << "\\\\ " << endl;
        (void)err;
	
	cout << "\\hline \\end{tabular} \\end{center} \\end{table} " << endl;
}
