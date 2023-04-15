#include<iostream>
#include<xmmintrin.h> //SSE
#include<emmintrin.h> //SSE2
#include<pmmintrin.h> //SSE3
#include<tmmintrin.h> //SSSE3
#include<smmintrin.h> //SSE4.1
#include<nmmintrin.h> //SSSE4.2
#include<immintrin.h> //AVX、AVX2
#include<windows.h>
#include<iomanip>
#include<ctime>
using namespace std;
void generateSample(float** A, int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			A[i][j] = 0;//下三角赋值为0;
		}
		A[i][i] = 1.0;//对角线赋值为1;
		for (int j = i; j < N; j++) {
			A[i][j] = rand();//上三角赋值为任意值;
		}
	}
	for (int k = 0; k < N; k++) {
		for (int i = k + 1; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i][j] += A[k][j];//每一行都加上比自己下标小的行;
			}
		}
	}
}
void show(float** A, int N) {//打印结果;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << fixed << setprecision(0) << A[i][j] << " ";
		}
		cout << endl;
	}
}

//串行算法:
void serialSolution(float** A, int N) {
	for (int k = 0; k < N; k++) {
		for (int j = k + 1; j < N; j++) {
			A[k][j] /= A[k][k];
		}
		A[k][k] = 1.0;
		for (int i = k + 1; i < N; i++) {
			for (int j = k + 1; j < N; j++) {
				A[i][j] -= A[i][k] * A[k][j];
			}
			A[i][k] = 0;
		}
	}
}

//并行算法:
void parallelSolution(float** A, int N) {
	for (int k = 0; k < N; k++) {
		__m128 vt = _mm_set1_ps(A[k][k]);
		for (int j = k + 1; j + 4 <= N; j += 4) {
			__m128 va = _mm_loadu_ps(&A[k][j]);
			va = _mm_div_ps(va, vt);
			_mm_storeu_ps(&A[k][j], va);
			if (j + 8 > N) {//处理末尾
				while (j < N) {
					A[k][j] /= A[k][k];
					j++;
				}
				break;
			}
		}
		A[k][k] = 1.0;
		for (int i = k + 1; i < N; i++) {
			__m128 vaik = _mm_loadu_ps(&A[i][k]);
			for (int j = k + 1; j + 4 <= N; j += 4) {
				__m128 vakj = _mm_loadu_ps(&A[k][j]);
				__m128 vaij = _mm_loadu_ps(&A[i][j]);
				__m128 vx = _mm_mul_ps(vakj, vaik);
				vaij = _mm_sub_ps(vaij, vx);
				_mm_storeu_ps(&A[i][j], vaij);
				if (j + 8 > N) {//处理末尾
					while (j < N) {
						A[i][j] -= A[i][k] * A[k][j];
						j++;
					}
					break;
				}
			}
			A[i][k] = 0;
		}
	}
}

int main() {
	float** A;
	float** B;
	int N = 1280;
	A = new float* [N];
	for (int i = 0; i < N; i++) {
		A[i] = new float[N];//申请空间;
	}
	B = new float* [N];
	for (int i = 0; i < N; i++) {
		B[i] = new float[N];//申请空间;
	}

	long long head_1, tail_1, freq;
	long long head_2, tail_2;
	int step = 64;//确定跨度;
	for (int n = step; n <= N; n += step) {//确定问题的规模

		generateSample(A, N);
		generateSample(B, N);
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

		//串行算法:
		QueryPerformanceCounter((LARGE_INTEGER*)&head_1);
		QueryPerformanceCounter((LARGE_INTEGER*)&tail_1);
		long long limit = 0.1 * freq;
		int counter_1 = 0;
		while ((tail_1 - head_1) < limit) {
			serialSolution(A, n);
			counter_1++;
			QueryPerformanceCounter((LARGE_INTEGER*)&tail_1);
		}

		//并行算法:
		QueryPerformanceCounter((LARGE_INTEGER*)&head_2);
		QueryPerformanceCounter((LARGE_INTEGER*)&tail_2);
		int counter_2 = 0;
		while ((tail_2 - head_2) < limit) {
			parallelSolution(B, n);
			counter_2++;
			QueryPerformanceCounter((LARGE_INTEGER*)&tail_2);
		}

		//时间统计:
		float seconds_1 = (tail_1 - head_1) * 1000.0 / freq;
		float seconds_2 = (tail_2 - head_2) * 1000.0 / freq;
		cout << fixed << setprecision(5);
		cout << "问题规模" << n << "	串行算法平均用时:" << seconds_1 / counter_1;
		cout << "	并行算法平均用时" << seconds_2 / counter_2 << endl;
	}
	return 0;
}