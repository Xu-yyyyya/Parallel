#include <iostream>
#include<windows.h>

using namespace std;

int arr[60000];
void chain_unroll(int arr[],int &sum,int N){
            for(int i =0;i<N;i++){
               sum+=arr[i];
            }
        sum=0;

}

void chain_two(int arr[],int &sum,int N){
            int sum1 = 0;
            int sum2 = 0;
            for (int i = 0;i < N; i += 2) {
                sum1 += arr[i];
                sum2 += arr[i + 1];
                }
            sum = sum1 + sum2;
            sum=0;
}
void recursion(int n){
    if(n==1)
        return;
    else{
        for(int i=0;i<n/2;i++)
            arr[i]+=arr[n-1-i];
        n=n/2;
        recursion(n);
    }
}

void xunhuan(int arr[],int N){
    for (int m = N; m > 1; m /= 2)
             for (int i = 0; i < m / 2; i++)
             arr[ i ] = arr[i*2] + arr[i*2 + 1];
}
int main()
{
    int N;
    cin>>N;
        int sum = 0;
        for(int i=0;i<N;i++){
            arr[i]=i;
        }
   /*//ƽ���㷨

        for(int num=0;num<4000;num++)
        chain_unroll(arr,sum,N);
        cout<<"ƽ���㷨�ĳ���ִ�����"<<endl;


     //����·ʽ
     for(int num=0;num<4000;num++)
      chain_two(arr,sum,N);

        cout<<"��·��ʽ�㷨ִ�����"<<endl;


    //β�ݹ�
            for(int num=0;num<4000;num++){
            recursion(N);
            for(int i=0;i<N;i++){
            arr[i]=i;
        }
            }

        cout<<"β�ݹ��㷨�ĳ���ִ�����"<<endl;
         for(int i=0;i<N;i++){
            arr[i]=i;
        }*/

        for(int num=0;num<4000;num++){
       xunhuan(arr,N);
       for(int i=0;i<N;i++){
            arr[i]=i;
        }
        }
        cout<<"����ѭ���㷨�ĳ���ִ�����"<<endl;

    return 0;
}
