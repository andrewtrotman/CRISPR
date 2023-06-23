//gather.c
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

#define N 1024
#define R 1000000

void foo_auto(double * restrict a, double * restrict b, int *idx, int n);
void foo_AVX2(double * restrict a, double * restrict b, int *idx, int n);
void foo_AVX512(double * restrict a, double * restrict b, int *idx, int n);
void foo1(double * restrict a, double * restrict b, int *idx, int n);
void foo2(double * restrict a, double * restrict b, int *idx, int n);
void foo3(double * restrict a, double * restrict b, int *idx, int n);


double test(int *idx, void (*fp)(double * restrict a, double * restrict b, int *idx, int n)) {
  double a[N];
  double b[N];
  double dtime;

  for(int i=0; i<N; i++) a[i] = 1.0*N;
  for(int i=0; i<N; i++) b[i] = 1.0;
  fp(a, b, idx, N);
  dtime = -omp_get_wtime();
  for(int i=0; i<R; i++) fp(a, b, idx, N);
  dtime += omp_get_wtime();
  return dtime;
}

int main(void) {

  //for(int i=0; i<N; i++) idx[i] = N - i - 1;
  //for(int i=0; i<N; i++) idx[i] = i;
  //for(int i=0; i<N; i++) idx[i] = rand()%N;

  //for(int i=0; i<R; i++) foo2(a, b, idx, N);
  int idx[N];
  double dtime;
  int ntests=2;
  void (*fp[4])(double * restrict a, double * restrict b, int *idx, int n);
  fp[0] = foo_auto;
  fp[1] = foo_AVX2;
#if defined ( __AVX512F__ ) || defined ( __AVX512__ )
  fp[2] = foo_AVX512;
  ntests=3;
#endif     

  for(int i=0; i<ntests; i++) { 
    for(int i=0; i<N; i++) idx[i] = 0;
    test(idx, fp[i]);
    dtime = test(idx, fp[i]);
    printf("%.2f      ", dtime);

    for(int i=0; i<N; i++) idx[i] = i;
    test(idx, fp[i]);
    dtime = test(idx, fp[i]);
    printf("%.2f      ", dtime);

    for(int i=0; i<N; i++) idx[i] = N-i-1;
    test(idx, fp[i]);
    dtime = test(idx, fp[i]);
    printf("%.2f      ", dtime);

    for(int i=0; i<N; i++) idx[i] = rand()%N;
    test(idx, fp[i]);
    dtime = test(idx, fp[i]);
    printf("%.2f\n", dtime);
  }

  for(int i=0; i<N; i++) idx[i] = 0;
  test(idx, foo1);
  dtime = test(idx, foo1);
  printf("%.2f      ", dtime);

  for(int i=0; i<N; i++) idx[i] = i;
  test(idx, foo2);
  dtime = test(idx, foo2);
  printf("%.2f      ", dtime);

  for(int i=0; i<N; i++) idx[i] = N-i-1;
  test(idx, foo3);
  dtime = test(idx, foo3);
  printf("%.2f      ", dtime);
  printf("NA\n");
}

//foo2.c
//#include <x86intrin.h>
#include <immintrin.h>
void foo_auto(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i++) b[i] = a[idx[i]];
}

void foo_AVX2(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i+=4) {
    __m128i vidx = _mm_loadu_si128((__m128i*)&idx[i]);
    __m256d av = _mm256_i32gather_pd(&a[i], vidx, 8);
    _mm256_storeu_pd(&b[i],av);
  }
}

#if defined ( __AVX512F__ ) || defined ( __AVX512__ )
void foo_AVX512(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i+=8) {
    __m256i vidx = _mm256_loadu_si256((__m256i*)&idx[i]);
    __m512d av = _mm512_i32gather_pd(vidx, &a[i], 8);
    _mm512_storeu_pd(&b[i],av);
  }
}
#endif

void foo1(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i++) b[i] = a[0];
}

void foo2(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i++) b[i] = a[i];
}

void foo3(double * restrict a, double * restrict b, int *idx, int n) {
  for(int i=0; i<n; i++) b[i] = a[n-i-1];
}
