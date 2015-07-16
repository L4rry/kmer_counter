#ifndef BP_BUS64_MPI_H_INCLUDED
#define BP_BUS64_MPI_H_INCLUDED

void workerCode(int id);

//////////////////////////////////////konfiguracne parametre
#define BUFSZ 1000000
// #define TARGET_PERIOD 1         // 1-millisecond target interrupt period
#define bitwisePartition
#define baseBitShift 8   // pocet bucketov load balancing 2^n  (BASEBITSHIFT je to uvodne delenie - je tam 8 -> 2exp8 = 256, odporucam 16*pocet jadier)
#define maxLgNBs 12				//max. pocet bucketov v deleni - 2^n (MAXLGNBSS sa tyka bucket sortu, da sa pskusat ci 11, 12, 13, 14 vyssie by som nesiel, ale nastavene je to na 12, ze to takdobre chodilo)

int NWS;        //pocet workerov
int NWThs;
int NMS;        //pocet masterov
int NMThs;
uint64_t DATALENG;

uint64_t *data;
uint64_t *aux;

typedef struct _callData {
    uint64_t *data;
    uint64_t *aux;
    int start;
    int count;
    int rMin;   //ak bitwisePart <= pocet bitov zhora
    int rMax;
    int *bCounts;
    int THBUFSZ;
    int NBUCKETS;
    pthread_mutex_t mutex;
} callData;




#endif // BP_BUS64_MPI_H_INCLUDED
