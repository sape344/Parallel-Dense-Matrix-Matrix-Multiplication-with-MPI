#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include <string>
#include <sys/time.h>
#include "mpi.h"
//#define SIZE 10000
#define Rand_low_bound 1.0
#define Rand_upper_bound 1000.0
typedef float deftype;


template<typename T>
T* createArray2DDynamic(const int size);
template<typename T>
T* createArray2DDynamic(const int n,const int m);
template<typename T>
void deleteArray2DDynamic(T* array, const int size);
template<typename T>
T* matMult2DDynamic( T* A,  T* B, const int size);
template<typename T>
T* matMult2DDynamic (T* A,  T* B, const int size,const int local_size) ;
template<typename T>
void fill2DRandomNumber(T* A,const int size );
template<typename T>
void Print2DDynamicArray( T* C, const int size);

template<typename T>
T* matMult1DDynamic(const T* A,const T* B, const int size_i,const int size_j);

void ArraysForScatterv(int size,int nproc,int*& sendcounts, int*& displs ){
    if(sendcounts==nullptr){
        sendcounts=new int[nproc];
    }        
    
    if(displs==nullptr){
        displs=new int[nproc];
    }

    for (int i = 0,sum=0,rem=(size%nproc)*size; i < nproc; i++)
    {
        sendcounts[i]=size*(size/nproc);
        if(rem>0 && nproc-1== i){
            sendcounts[i]+=rem;
            rem=0;
        }
        displs[i] = sum;
        sum += sendcounts[i];
    }
    

}
    int nproc,rank;


int main(int argc, char **argv){
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int SIZE=std::atoi(argv[1]);
    //const int SIZE=25000;
    
    int local_size=SIZE/nproc;

    


    deftype *A=nullptr;
    deftype *A_local=nullptr;
    deftype *B_local=nullptr;
    deftype *B=nullptr;
    deftype *C=nullptr;
    int rem= SIZE%nproc;
    int *sendcounts=nullptr;
    int *displs=nullptr;



    

    ArraysForScatterv(SIZE,nproc,sendcounts, displs );
    
    double t1,t2;
    if(rank==nproc-1){
        local_size=SIZE-local_size*(nproc-1);
    }    
    B = createArray2DDynamic<deftype>(SIZE);

    if (rank==0)
    {       
        C = createArray2DDynamic<deftype>(SIZE);
        A = createArray2DDynamic<deftype>(SIZE);        
        fill2DRandomNumber(A,SIZE);
        fill2DRandomNumber(B,SIZE);
      //  Print2DDynamicArray(B,SIZE);

        
    }
  
    if (rank==0)
    {
        t1 = MPI_Wtime(); 
    }


    
    A_local= createArray2DDynamic<deftype>(SIZE,local_size);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(A,sendcounts,displs,MPI_FLOAT,A_local,local_size*SIZE,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(B,SIZE*SIZE,MPI_FLOAT,0,MPI_COMM_WORLD);
    auto C_local= matMult2DDynamic( A_local, B ,local_size, SIZE);
    MPI_Gatherv(C_local,SIZE*local_size,MPI_FLOAT,C,sendcounts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
    {
        t2 = MPI_Wtime(); 

        std::cout<<"Core Number= "<<nproc<<" Matrix Size= "<<SIZE<< " Wall time= "<<t2-t1<<"\n";

    }
    
    //Print2DDynamicArray(C,SIZE);
    
    deleteArray2DDynamic(B, SIZE);
    deleteArray2DDynamic(A_local, local_size);
    deleteArray2DDynamic(C_local, local_size);

    if (rank==0){  
        deleteArray2DDynamic(A, SIZE);
        deleteArray2DDynamic(C, SIZE);
     }


    MPI_Finalize();

            
    return 0;
}






template<typename T>
T* matMult2DDynamic( T* A,  T* B, const int size_i,const int size_j) {
    auto C = createArray2DDynamic<deftype>(size_i,size_j);
    deftype sum;
    deftype r;


    for ( int i = 0; i < size_i*size_j; i++)
    {
        C[i]=0;
    }
    
        for (int i = 0; i < size_i; i++)
            {
            for (int k = 0; k < size_j; k++)     
                {
                    r = A[i*size_j+k];

                     for (int j = 0;j < size_j; j++)
                    {
                        C[i*size_j+j] += B[k*size_j+j] * r;

                    }

                }
            }


    return C;
}

template<typename T>
T* matMult1DDynamic(const T* A,const T* B, const int size_i,const int size_j){
     auto C = createArray2DDynamic<deftype>(size_i,size_j);

    for(int i=0; i<size_i;i++){
        for(int j=0; j<size_j;j++){
            C[i*size_i+j]=0;
            for (int k = 0; k < size_j; k++)
            {
                //C[size_i*i+j]+=A[i*size_i+k]*B[j+size_j*k];
                C[i*size_i+j]+=A[i*size_i+k]*B[j+size_j*k];
            }
        } 
    } 
    

    return C;

}

template<typename T>
T* createArray2DDynamic(const int size) {
    T* AAA;
    AAA = new T [size*size];

    return AAA;
}

template<typename T>
T* createArray2DDynamic(const int m,const int n) {
    T* AAA;
    AAA = new T [n*m];

    return AAA;
}


template<typename T>
void Print2DDynamicArray( T* C, const int size){
    std::cout<<"\n--------------------------------------"<<std::endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {       

                std::cout<<C[i*size+j]<<"\t";
            
            }
                std::cout<<"\n";

        }
        std::cout<<"\n--------------------------------------\n";
}



template<typename T>
void fill2DRandomNumber(T* A,const int size ){
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<T> dist6(Rand_low_bound,Rand_upper_bound);
    //std::uniform_int_distribution<T> dist6(Rand_low_bound,Rand_upper_bound);
   

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i*size+j]=dist6(rng);
        }
        
    }


}

template<typename T>
void deleteArray2DDynamic(T* array, const int size) {
if(array != nullptr){
    delete[] array;

}
}

