#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <inttypes.h>
#include <pthread.h>
#include <semaphore.h>
#include <mpi.h>
#include <string.h>

#include "BP_BuS64_MPI.h"

typedef struct _MScallData {
    int start;
    int count;
    sem_t *semaf;
} MScallData;

typedef struct _keyVal {
  uint64_t kmer;
  int count;
} keyVal;


short myID;
int myBaseBitShift;


void prove(uint64_t *arr, int count) {
    int i;
    for(i=1; i<count; i++) {
        if(arr[i-1]>arr[i]){
            fprintf(stderr, "chyba %d %d %llu %llu\n", myID, i, arr[i-1], arr[i]);
            return;
        }
    }
}

void insertion_sort(uint64_t *a, int n) {
    int i, j;
    uint64_t value;
    for (i = 1; i < n; i++) {
        value = a[i];
        for (j = i; j > 0 && value < a[j - 1]; j--) {
            a[j] = a[j - 1];
        }
        a[j] = value;
    }
}

void quick_sort (uint64_t *data, int *pivInxs, int n) {
    int i, j, k;
    uint64_t pivot, auxn;
    int npivsL = 0;
    int npivsR = 0;

    if (n < 32){
        insertion_sort(data, n);
        return;
    }

    pivot = data[n / 2];
    i = 0;
    j = n - 1;
    for (; ; i++, j--) {
        while (data[i] < pivot){
            i++;
        }
        while (pivot < data[j]){
            j--;
        }
        if (i >= j)
            break;
       
        if(data[j]==pivot){
            pivInxs[npivsL] = i;
            npivsL++;
        }
        if(data[i]==pivot){
            npivsR++;
            pivInxs[n-npivsR] = j;
        }

        auxn = data[i];
        data[i] = data[j];
        data[j] = auxn;
    }

    int midInx = i-1;
    for(npivsL--;npivsL>=0;npivsL--){
        if(data[midInx]!=pivot){
            auxn = data[midInx];
            data[midInx] = data[pivInxs[npivsL]];
            data[pivInxs[npivsL]] = auxn;
        }
        midInx--;
    }
    int leftEnd = midInx+1;


    midInx = i;
    for(npivsR--;npivsR>=0;npivsR--){
        if(data[midInx]!=pivot){
            auxn = data[midInx];
            data[midInx] = data[pivInxs[n-npivsR-1]];
            data[pivInxs[n-npivsR-1]] = auxn;
        }
        midInx++;
    }
    int rightStart = midInx;

    quick_sort(data, pivInxs, leftEnd);
    quick_sort(data + rightStart, pivInxs + rightStart, n - rightStart);
}


void bucketSortSerial(uint64_t *data, uint64_t *aux, int dataCount, int _baseBitShift, int lgNBUCKETS, short beNested) { //vstup v aux, vystup v data   // dataCount - pocet prvkov, _bbs - pocet zhora rovnakych bitov, lgNB - log. pocetu vedierok

    int i;
    
    int NBUCKETS = (1 << lgNBUCKETS);
    int bitShift = 64 - lgNBUCKETS;

    int bCounts[NBUCKETS];       memset(bCounts, 0, NBUCKETS*sizeof(int));
    int cumBCounts[NBUCKETS];    memset(cumBCounts, 0, NBUCKETS*sizeof(int));//aka bucket start

   
    ////////////////////////////////////////////////////1st pass (zisti pocetnosti v aux)
    {
        for(i=0; i<dataCount; i++) {
            int inx = (aux[i] << _baseBitShift) >> (bitShift);
            bCounts[inx]++;
        }
        cumBCounts[0] = 0;
        for(i=1; i<NBUCKETS; i++) {
            cumBCounts[i] = cumBCounts[i-1] + bCounts[i-1] ;
        }
    }

    ////////////////////////////////////////////////////2nd pass (prekopiruj z aux na prislusne miesto v data)
    {
        for(i=0; i<dataCount; i++) {
            int bIdx = (aux[i] << _baseBitShift) >> (bitShift);
            data[cumBCounts[bIdx]] = aux[i];
            cumBCounts[bIdx]++;
        }
        for(i=NBUCKETS-1; i>0; i--) {
            cumBCounts[i] = cumBCounts[i-1];
        }
        cumBCounts[0] = 0;
    }


    ////////////////////////////////////////////////////bucket sort
    {
        if(beNested == 0){
            for(i=0; i<NBUCKETS; i++) {
               quick_sort(&data[cumBCounts[i]], (int*)&aux[cumBCounts[i]], bCounts[i]);
            }
        }else{
            for(i=0; i<NBUCKETS; i++) {
              int n = floor(log(bCounts[i] / 25) / log(2) );
              int nextShift = n < maxLgNBs ? n:maxLgNBs;  /// parameter - max. pocet bucketov v 2. deleni - 2^n
              
              if(nextShift<4){
                quick_sort(&data[cumBCounts[i]], (int*)&aux[cumBCounts[i]], bCounts[i]);
                memcpy(&aux[cumBCounts[i]], &data[cumBCounts[i]], 8*bCounts[i]); //lebo vysledok cakam v globalnom data, co je teraz aux
              }
              else{
                bucketSortSerial(&aux[cumBCounts[i]], &data[cumBCounts[i]], bCounts[i], _baseBitShift+64-bitShift, nextShift, 0);
              }
            }
        }
    }
}

void workerCode(int id) { 
		printf("WORKER START\n");
    struct timeval start, end;
    double timeSpent;
    myID = id;
    uint64_t *data[NWThs];
    uint64_t *aux[NWThs];
    int i,j;

    MPI_Status recvStatus;
    int dataSize[NWThs]; //celkova velkost pre kazdy thread
    int dataLevelCount[NWThs];   //pocet hladin pre kazdy thread
    int *dataLevelSizes[NWThs];   //velkost hladin pkt
    int *cumDataLevelSizes[NWThs];      //pociatocny index pkt
    
    ///////////////////////////////////////////////////////////inicializacia
    int mesAux[2 * NWThs];   //(pocet, velkost, pocet, velkost, ...)
    MPI_Recv((void*)mesAux, 2 * NWThs, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &recvStatus);
    printf("Worker %d - TAG 2 passed\n", id);
//     printf("<<<\n");
//     for (i = 0; i < 2 * NWThs; i++) {
//     	printf("%d\n", mesAux[i]);
//     }
//     printf("<<<\n");
    
    for(i = 0;i < NWThs; i++){          //nastav velkosti a pocty hladin
        dataLevelCount[i] = mesAux[i*2];
        dataSize[i] = mesAux[i*2+1];
				printf("dataLevelCount[%d] = %d\n", i, mesAux[i*2]);
        printf("dataSize[%d] = %d\n", i, mesAux[i*2+1]);
        dataLevelSizes[i] = (int*) malloc(dataLevelCount[i]*sizeof(int));
        cumDataLevelSizes[i] = (int*) malloc(dataLevelCount[i]*sizeof(int));
    }

    //prijmi velkosti hladin
    printf("Worker %d - TAG 3 started\n", id);
    for (i = 0; i < NWThs; i++){
        MPI_Recv((void*)dataLevelSizes[i], dataLevelCount[i], MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &recvStatus);
    }
    printf("Worker %d - TAG 3 passed\n", id);
        
    //nastav start indexy
    for (i = 0; i < NWThs; i++) {
        cumDataLevelSizes[i][0] = 0;
        
        for(j = 1; j < dataLevelCount[i]; j++){
          cumDataLevelSizes[i][j] = cumDataLevelSizes[i][j-1] + dataLevelSizes[i][j-1];
        }
    }       
    //alokuj data
    for (j = 0;j < NWThs; j++) {
        data[j] = (uint64_t*)malloc(dataSize[j] * 8);
    } 
// 
//  		int *hlp;
//     MPI_Send(hlp, 1, MPI_INT, 0, 444, MPI_COMM_WORLD);
// 
    ///////////////////////////////////////////////////////prijem dat
		printf("Worker %d - TAG 0 & 1 started\n", id);
    int endedTmss = 0;
    while(1) {

        int dataLMes[3];
        MPI_Recv((void*)dataLMes, 3, MPI_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvStatus);
        int size = dataLMes[0];
        int thread = dataLMes[1];
        int inx = dataLMes[2];
        if(size==0) {
            endedTmss++;
            if(endedTmss == NMS*NWThs) {
                break;
            }
            continue;
        }
        uint64_t *destData = data[thread];  
        MPI_Recv((void*)&(destData[cumDataLevelSizes[thread][inx]]), size, MPI_UNSIGNED_LONG_LONG, recvStatus.MPI_SOURCE, 1, MPI_COMM_WORLD, &recvStatus );
        cumDataLevelSizes[thread][inx] += size;
    }
    printf("Worker %d - TAG 0 & 1 passed\n", id);
   // MPI_Barrier(MPI_COMM_WORLD);
    
    //nainicializuj naspat start indexy
    for(i = 0; i < NWThs; i++) {
      for(j = dataLevelCount[i]-1; j > 0; j--) {
        cumDataLevelSizes[i][j] = cumDataLevelSizes[i][j-1];
      }
      cumDataLevelSizes[i][0] = 0;
    }

    //zalokuj pomocne pole
    for(j = 0; j < NWThs; j++){
        aux[j] = (uint64_t*)malloc(dataSize[j]*8);
    }
  
    gettimeofday(&start, NULL);

    /// /////////////////////////////////////////////////////pocitanie
    
    #pragma omp parallel num_threads(NWThs)
    {
        int thInx = omp_get_thread_num();
        int k;
        uint64_t *myData = data[thInx];
        uint64_t *myAux = aux[thInx];    
        for(k=0;k<dataLevelCount[thInx];k++){
            bucketSortSerial(&(myAux[cumDataLevelSizes[thInx][k]]), &(myData[cumDataLevelSizes[thInx][k]]), dataLevelSizes[thInx][k], baseBitShift, maxLgNBs, 1); 
        }               
    }

    gettimeofday(&end, NULL);
    timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;

    prove(data[0], dataSize[0]);
    prove(data[1], dataSize[1]);
    
    // ZAPIS SORTED
//     for(j=0;j<NWThs;j++){
// //         printf("zapis sorted: ...\n");
//         char oFileName[200];
// //         sprintf(oFileName, "/work/projects/3DNA-2013/data/freqCountData/n8/res/IT_Met_1a.chr22.sorted.part%d.bin", (id-NMS)*NWThs+j);
// 				sprintf(oFileName, "/work/projects/3DNA-2013/data/freqCountData/data/SRX040485/sorted/bigger/aa.part%d.bin", (id-NMS)*NWThs+j);
//         FILE *fout = fopen(oFileName, "wb");
//         if(fout==NULL){
//           printf("chyba otvarania w suboru %d %d\n", id, j);
//           exit(1);
//         }
//         int wrote = fwrite(data[j], 8, dataSize[j], fout);
//         fflush (fout);
//         printf("%d z %d zapisane do %s\n", wrote, dataSize[j], oFileName);
//         fclose(fout);
//     }
    
    int length;
    char hostID[100];
    MPI_Get_processor_name(hostID, &length);  
    for(j = 0; j < NWThs; j++){
        printf("worker %d na %s thread %d spracoval %d od %llu do %llu za %lf\n", id, hostID, j, dataSize[j], data[j][0], data[j][dataSize[j]-1], timeSpent          );       
    } 
    
    #pragma omp parallel num_threads(NWThs)
    {
        struct timeval start, end;
        
        gettimeofday(&start, NULL);
        int thInx = omp_get_thread_num();
        keyVal *myRes = (keyVal*) aux[thInx];
        uint64_t *myData = data[thInx];
        
        int i;
        
        uint64_t curKmer = myData[0];
        int curCount = 1;
        int distCount = 0;
        
        for(i=1; i<dataSize[thInx];i++){
          if(myData[i] != curKmer){
            myRes[distCount].kmer = curKmer;
            myRes[distCount].count = curCount;
            distCount++;
            curKmer = myData[i];
            curCount=1;          
          }else{
            curCount++;
          }        
        } 
        gettimeofday(&end, NULL);
        double timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;    
        
        
        printf("counting: %d %d %d %lf\n", id, thInx, distCount, timeSpent);
//         
//         // ZAPIS COUNTED
// //         printf("zapis counted...\n");
//         char oFileName[200];
// //         sprintf(oFileName, "/work/projects/3DNA-2013/data/freqCountData/n8/res/IT_Met_1a.chr22.counted.part%d.bin", id*NWThs + thInx);
// 				sprintf(oFileName, "/work/projects/3DNA-2013/data/freqCountData/data/SRX040485/counted/bigger/aa.counted.part%d.bin", id*NWThs + thInx);
//         FILE *fout = fopen(oFileName, "wb");
//         if(fout==NULL){
//           printf("chyba otvarania w suboru %d %d\n", id, j);
//           exit(1);
//         }
//         int wrote = fwrite(myRes, 12, distCount, fout);
//         fflush (fout);
//         fclose(fout);  
        
    }

    
    
}