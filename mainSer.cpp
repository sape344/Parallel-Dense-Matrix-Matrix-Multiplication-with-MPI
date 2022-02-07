#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include <string>
#include <sys/time.h>
///#include "mulmat.h"
#define N 1000
#define Rand_low_bound 1.0
#define Rand_upper_bound 1000.0


typedef float deftype; // defined worked type 
enum indexRearraged { ijk, jik, jki, kji, kij, ikj };

template<typename T>
T** createArray2DDynamic(const int size);
template<typename T>
void deleteArray2DDynamic(T** array, const int size);
template<typename T>
T** matMult2DDynamic( T** A,  T** B, const int size);


template<typename T>
T** matBlockMult2DDynamic( T** A,  T** B, const int size, const int blockSize);
template<typename T>
T** matMult2DDynamic( T** A,  T** B, const int size,const indexRearraged type_loop);
template<typename T>
void fill2DRandomNumber(T** A,const int size );
template<typename T>
void Print2DDynamicArray( T** C, const int size);




int main(int argc, char **argv){
    
    int size, block_size;
    indexRearraged type;
  
    if (argc==2){
        size= std::stoi(argv[1]);
        //type=static_cast<indexRearraged>(std::stoi(argv[2]));
        //block_size=100;
    }


    auto A = createArray2DDynamic<deftype>(size);
    auto B = createArray2DDynamic<deftype>(size);
    fill2DRandomNumber(A,size);
    fill2DRandomNumber(B,size);
    
    auto C1= matMult2DDynamic<deftype>(A,B,size,static_cast<indexRearraged>(5));    
    //Print2DDynamicArray<deftype>(A,size);
    deleteArray2DDynamic(A, size);
    deleteArray2DDynamic(B, size);
    deleteArray2DDynamic(C1, size);

            
    return 0;
}




template<typename T>
T** matMult2DDynamic( T** A,  T** B, const int size,const indexRearraged type_loop) {
    auto C = createArray2DDynamic<deftype>(size);
    deftype sum;
    deftype r;

    struct timeval begin, end;    
 
    switch (type_loop)
    {
    case ijk:
        gettimeofday(&begin, 0);
        for (int i = 0; i < size; i++)
            {
            for (int j = 0; j < size; j++)
                {
                    sum= 0;
                    //#pragma clang loop unroll(enable)
                    for (int k = 0; k < size; k++)
                    {
                        sum += A[i][k] * B[k][j];
                    }
                    C[i][j]=sum;
                }
            }
        gettimeofday(&end, 0);                
        break;

    case jik:
        gettimeofday(&begin, 0);     
        for (int j = 0; j < size; j++)        
            {
            for (int i = 0; i < size; i++)
                {
                    sum = 0;
                    //#pragma clang loop unroll(enable)
                    for (int k = 0; k < size; k++)
                    {
                        sum += A[i][k] * B[k][j];
                    }
                    C[i][j]=sum;
                }
            }
        gettimeofday(&end, 0);        
        break;
    case jki:
        gettimeofday(&begin, 0);
        for (int j = 0; j < size; j++)        
            {
            for (int k = 0; k < size; k++)
                {
                    r = B[k][j];
                    //#pragma clang loop unroll(enable)
                     for (int i = 0; i < size; i++)
                    {
                        C[i][j] += A[i][k] * r;
                    }

                }
            }        
        gettimeofday(&end, 0);
        break;
    case kji:       
        gettimeofday(&begin, 0);
        for (int k = 0; k < size; k++)        
            {
            for (int j = 0; j < size; j++)
                {
                    r = B[k][j];
                    for (int i = 0; i < size; i++)
                    {
                        C[i][j] += A[i][k] * r;
                    }

                }
            }
        gettimeofday(&end, 0);        
        break;
    case kij:
        gettimeofday(&begin, 0);
        for (int k = 0; k < size; k++)        
            {
            for (int i = 0; i < size; i++)
                {
                    r = A[i][k];
                     for (int j = 0;j < size; j++)
                    {
                        C[i][j] += B[k][j] * r;
                    }

                }
            }        
        gettimeofday(&end, 0);
        break;

    case ikj:
        gettimeofday(&begin, 0);
        for (int i = 0; i < size; i++)
            {
            for (int k = 0; k < size; k++)     
                {
                    r = A[i][k];
                    //#pragma clang loop unroll(enable)
                     for (int j = 0;j < size; j++)
                    {
                        C[i][j] += B[k][j] * r;
                    }

                }
            }
        gettimeofday(&end, 0);        
        break;
    
    default:
        gettimeofday(&begin, 0);
        gettimeofday(&end, 0);    

        break;
    }
    std::cout<<"Simple Matrix multiplication WallTime(s): " <<(end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec)*1e-6)<<" Size: "<<size << " index seq: "<<int(type_loop)<< std::endl;



    return C;
}
template<typename T>
T** matBlockMult2DDynamic( T** A,  T** B, const int size, const int blockSize) {
    auto C = createArray2DDynamic<deftype>(size);
    struct timeval begin, end;
    int sum;
    gettimeofday(&begin, 0);
    for(int i=0;i<size;i++){
        //#pragma clang loop unroll(enable)
        for(int j=0;j<size;j++){
            C[i][j]=0.0;        
        } 
    }

    for(int bi=0; bi<size; bi+=blockSize){
        for(int bk=0; bk<size; bk+=blockSize){
            for(int bj=0; bj<size; bj+=blockSize){
                for(int i=0; i<  ((bi+blockSize) >size ? size-bi:(blockSize)); i++){
                    for(int k=0; k<  ((bk+blockSize) > size?size-bk:(blockSize)); k++){
                        //#pragma clang loop unroll(enable)
                        for(int j=0; j<  ((bj+blockSize) > size ? size-bj:(blockSize)); j++){
                            C[bi+i][bj+j] += A[bi+i][bk+k]*B[bk+k][bj+j];
                        }
                    }
                }
            }
        }
    }



    gettimeofday(&end, 0);
    std::cout<<"Block Matrix multiplication WallTime(ms): " <<(end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec)*1e-6)*1e3<<" Size: "<<size << " Block Matrix size: "<<blockSize<< std::endl;

    return C;
}




template<typename T>
T** createArray2DDynamic(const int size) {
    T** AAA;
    AAA = new T* [size];
    for (int i = 0; i < size; i++)
    {
        AAA[i] = new T[size];
    }

    return AAA;
}

template<typename T>
void Print2DDynamicArray( T** C, const int size){
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {       

                std::cout<<C[i][j]<<"\t";
            
            }
                std::cout<<"\n";

        }
}



template<typename T>
void fill2DRandomNumber(T** A,const int size ){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<T> dist6(Rand_low_bound,Rand_upper_bound);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j]=dist6(rng);
        }
        
    }


}

template<typename T>
void deleteArray2DDynamic(T** array, const int size) {
    for (int i = 0; i < size; i++) {
        delete[] array[i];
    }
    delete[] array;
}

template<typename T>
T** matMult2DDynamic( T** A,  T** B, const int size) {
    auto C = createArray2DDynamic<deftype>(size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < size; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


