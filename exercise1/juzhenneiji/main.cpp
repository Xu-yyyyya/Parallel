#include <iostream>
#include<windows.h>

using namespace std;


float arr[10000][10000];

float sum[10000];

float a[10000];

int main()
{
        int N;
    while(cin>>N){


        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++)
                arr[i][j] =i*1.1+j*1.1;
            a[i]=1.1;
        }


        //平凡算法

 /*       long long head,tail,freq;
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        for(int num=0;num<10;num++){
            for(int i =0;i<N;i++){
                sum[i]=0.0;
                for(int j=0;j<N;j++)
                    sum[i] += a[i]*arr[j][i];
            }
        }

        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout<<"当问题规模N="<<N<<"时，平凡算法的程序执行时间为   "<<(tail-head)*100.0/freq<<"ms"<<endl;*/

        //Cache优化
        long long head1,tail1,freq1;
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        for(int num=0;num<10;num++){
            for( int i = 0; i < N; i++)
                 sum[i] = 0.0;
                 for(int j = 0; j < N; j++)
                 for(int i = 0; i < N; i++)
                 sum[i] += arr[j][i]*a[j];
        }

        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        cout<<"当问题规模N="<<N<<"时，cache优化算法的程序执行时间为   "<<(tail1-head1)*100.0/freq1<<"ms"<<endl;

    }
    return 0;
}
