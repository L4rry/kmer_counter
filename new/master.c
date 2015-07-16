#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
// #include <inttypes.h>
#include <string.h>
#include <unistd.h> 

#include "BP_BuS64_MPI.h"

#define BASEBITSHIFT 8
#define MAX_MEMORY_PROCESSOR_FILE_LOAD 4*134217728 //1073741824
#define RESERVE 1024
#define PART_SIZE_MAX 1073741824
#define KMER_SIZE 8 // 8B

MPI_Comm Comm_masters;
char BASES[255], BASES_2[255];

typedef struct FileReadStartEnd {
    size_t start;
    size_t end;
} FileReadStartEnd;

typedef struct PartLimit {
    int id;
    int start;
    int end;
    int num_parts;
    uint64_t size;
    uint64_t count;
} PartLimit;

void fillArray(uint64_t *data, int count, int id) { // count == DATALENG
	char oFileName[200];
	sprintf(oFileName, "/Users/larry/Desktop/tests/aa/aa_part%d.bin", id);
	FILE *fin = fopen(oFileName, "rb");

	if (fin == NULL) {
		exit(1);
	}
  fread(data, 8, count, fin);
  fclose(fin);
  printf("%s\n", oFileName);    
}

size_t getEndOfSequenceForward(int INPUT_FILE_TYPE, char *buf, size_t size_to) { // 0: raw, 1: fasta
    size_t i;
    size_t end_offset = 0;
    for (i = 0; i < size_to; i++) {
        end_offset++;
        if (buf[i] == '\n') {
            return end_offset;
        }
    }
    return 0;
}

void getPartsStartEnd(FILE *f, size_t f_size, int parts_count, size_t part_size, FileReadStartEnd *fileReadsStartsEnds) {
    size_t start = 0, end = 0, offset = 0;
    int i;
    char *buf = (char*)malloc(sizeof(char) * RESERVE);
    for (i = 0; i < parts_count; i++) {
        start = i * part_size + offset; // + 1 ? (okrem i == 0) ?
        end = (i+1) * part_size;
        fseek(f, (i+1) * part_size, SEEK_SET);
        fread(buf, sizeof(char), RESERVE, f);
        offset = getEndOfSequenceForward(1, buf, RESERVE);
        end += offset;
        memset(buf, '\0', RESERVE);
        if (end > f_size) {
            end = f_size;
        }
        fileReadsStartsEnds[i].start = start;
        fileReadsStartsEnds[i].end   = end;
    }
}

void getFileReadStartsEnds(int NMS, int myrank, char *input_filename, FileReadStartEnd **fileReadsStartsEnds, int *parts_count, size_t *part_size) {
//     printf("NMS = %d, myrank = %d, input_filename = %s\n", NMS, myrank, input_filename);
    // najst konce riadkov, a posunut starty/konce pre citanie zo suboru
    MPI_Status m_status;

    FILE *f = fopen(input_filename, "r");
    fseek(f, 0, SEEK_END);
    size_t f_size = ftell(f);
    fseek(f, 0, SEEK_SET);
    
    (*part_size)   = MAX_MEMORY_PROCESSOR_FILE_LOAD - RESERVE;
    (*parts_count) = f_size % (*part_size) == 0 ? f_size / (*part_size)  : f_size / (*part_size)  + 1;

    if ((*parts_count) <= NMS) {
        (*part_size) = f_size / NMS;
        (*parts_count) = NMS;
    }
    (*fileReadsStartsEnds) = (FileReadStartEnd*)malloc(sizeof(FileReadStartEnd) * (*parts_count));

    int i;
    size_t *_parts_read_start = (size_t*)malloc(sizeof(size_t) * (*parts_count));
    size_t *_parts_read_end   = (size_t*)malloc(sizeof(size_t) * (*parts_count));
    size_t *parts_read_start  = (size_t*)malloc(sizeof(size_t) * (*parts_count));
    size_t *parts_read_end    = (size_t*)malloc(sizeof(size_t) * (*parts_count));
    for (i = 0; i < (*parts_count); i++) {
        _parts_read_start[i] = 0;
        _parts_read_end[i] = 0;
        parts_read_start[i] = 0;
        parts_read_end[i] = 0;
    }

    // odoslat ostatnym informaciu o koncoch/zaciatkoch
    if (myrank == 0) {
        getPartsStartEnd(f, f_size,  (*parts_count), (*part_size), (*fileReadsStartsEnds));

        for (i = 0; i < (*parts_count); i++) {
            _parts_read_start[i] = (*fileReadsStartsEnds)[i].start;
            _parts_read_end[i] = (*fileReadsStartsEnds)[i].end;
        }

//        for (i = 0; i < (*parts_count); i++) {
//            printf("%d: start = %lu, end = %lu\n", i, (*fileReadsStartsEnds)[i].start, (*fileReadsStartsEnds)[i].end);
//        }

        if (NMS > 1) {
            for (i = 1; i < NMS; i++) {
                printf("[%d]: send\n", myrank);
                MPI_Send((void*)&parts_count[0], 1, MPI_INT, (int)i, 0, Comm_masters);
            }
        }
    } else { // prijat
        
        if (NMS > 1) {
                printf("[%d]: recv\n", myrank);
            MPI_Recv((void*)&parts_count[0], 1, MPI_INT, 0, 0, Comm_masters, &m_status);
        }
        
    }
//     MPI_Barrier(Comm_masters);
// 		printf("[%d]: parts_count = %d\n", myrank, (*parts_count));
// 		MPI_Finalize();
// 		exit(0);

    // spojit
    MPI_Allreduce(_parts_read_start, parts_read_start, (*parts_count), MPI_UINT64_T, MPI_SUM, Comm_masters);
    MPI_Allreduce(_parts_read_end, parts_read_end, (*parts_count), MPI_UINT64_T, MPI_SUM, Comm_masters);

    // ulozit spojene informacie o startoch/koncoch citania zo suboru
    if (myrank != 0) {
        for(i = 0; i < (*parts_count); i++) {
            (*fileReadsStartsEnds)[i].start = parts_read_start[i];
            (*fileReadsStartsEnds)[i].end = parts_read_end[i];
        }
    }
    free(_parts_read_start);
    free(_parts_read_end);
    free(parts_read_start);
    free(parts_read_end);
 
    fclose(f);
}
// z bufferu (celeho dlheho nacitaneho stringu zo suboru) spravi k-tice a spocitat aj n-tice
void doKmersNmers(char *buf, size_t buf_size, uint64_t *kmers, int *kmers_count, int KMER_LENGTH, int *nmers, int nBits) {
    size_t i;
    uint64_t kmer = 0;
    int kmer_length_counter = 0;
    int bitShift = 64 - nBits;

    for (i = 0; i < buf_size; i++) {

        if (buf[i] != '\n' && buf[i] != '\0' && buf[i] != EOF && buf[i] != 'N') {
            
            if (kmer_length_counter < KMER_LENGTH) {
                kmer |= (((uint64_t)BASES[buf[i]]) << (62 - kmer_length_counter * 2));
                kmer_length_counter++;
                if (kmer_length_counter == KMER_LENGTH) {
                    kmers[(*kmers_count)] = kmer;
                    nmers[kmer >> bitShift]++;
                    (*kmers_count)++;
                    //                    printf("first kmer = %llu\n", kmer);
                }
            } else {
                kmer = (uint64_t)((kmer << 2) | ((uint64_t)BASES_2[buf[i]])); // << 2
                kmers[(*kmers_count)] = kmer;
                nmers[kmer >> bitShift]++;
                (*kmers_count)++;
                //                printf("kmer = %llu\n", kmer);
            }
        } else {
            kmer = 0;
            kmer_length_counter = 0;
        }
    }
    //    printf("kmers_count = %d\n", *kmers_count);
}

void getKmersStartEnd(char *buf, size_t buf_size, uint64_t *kmers, int KMER_LENGTH, int start, int end, int nBits, int *kmers_count) {
    size_t i;
    uint64_t kmer = 0;
    int kmer_length_counter = 0;
    int bitShift = 64 - nBits;

    for (i = 0; i < buf_size; i++) {

        if (buf[i] != '\n' && buf[i] != '\0' && buf[i] != EOF && buf[i] != 'N') {
            
            if (kmer_length_counter < KMER_LENGTH) {
                kmer |= (((uint64_t)BASES[buf[i]]) << (62 - kmer_length_counter * 2));
                kmer_length_counter++;
                
                if (kmer_length_counter == KMER_LENGTH) {
                
                	if ( (kmer >> bitShift) >= start && (kmer >> bitShift) <= end ) {
                    kmers[(*kmers_count)] = kmer;
                    (*kmers_count)++;
                  }
//                    printf("first kmer = %llu\n", kmer);
                }
            } else {
                kmer = (uint64_t)((kmer << 2) | ((uint64_t)BASES_2[buf[i]])); // << 2     //TODO: vsade -> mensia k-tica ako 31 tu bude zle

                if ( (kmer >> bitShift) >= start && (kmer >> bitShift) <= end ) {
	                kmers[(*kmers_count)] = kmer;
  	              (*kmers_count)++;
  	            }
//                printf("kmer = %llu\n", kmer);
            }
        } else {
            kmer = 0;
            kmer_length_counter = 0;
        }
    }
    //    printf("kmers_count = %d\n", *kmers_count);
}


void getKmers(MPI_File m_f, FileReadStartEnd *fileReadsStartsEnds, int part_id, int start_level, int end_level, uint64_t *kmers, int KMER_LENGTH, int nBits, int *kmers_count) {
	 	printf("file read from %lu to %lu\n", fileReadsStartsEnds[part_id].start, fileReadsStartsEnds[part_id].end);
		int bytes_to_read = (int)(fileReadsStartsEnds[part_id].end - fileReadsStartsEnds[part_id].start);
		int bytes_really_read = 0;
		char *buf = (char*)malloc(sizeof(char) * bytes_to_read);

		MPI_Status m_status;
		MPI_File_seek(m_f, fileReadsStartsEnds[part_id].start, MPI_SEEK_SET);
		MPI_File_read(m_f, buf, bytes_to_read, MPI_CHAR, &m_status);
		MPI_Get_count(&m_status, MPI_CHAR, &bytes_really_read);

// 		void getKmersStartEnd(char *buf, size_t buf_size, uint64_t *kmers, int KMER_LENGTH, int start, int end, int nBits) {
		getKmersStartEnd(buf, bytes_really_read, kmers, KMER_LENGTH, start_level, end_level, nBits, kmers_count);
		free(buf);
}


void masterCode(int id) {
    int i, j;
    struct timeval start, end, preStart, wtStart, wtEnd;
    double timeSpent, time_spent, timeHist;
    
    int nBases = 4;
    int nBits = nBases * 2;
    int nParts = (1 << BASEBITSHIFT);
    
    int NBUCKETS = (1 << baseBitShift); // 8 => 256 (nParts)
    int bitShift = 64 - baseBitShift;   // 8 => 56

    data = (uint64_t*)malloc(DATALENG * sizeof(uint64_t));
    if (data == NULL) { fprintf(stderr, "NEED MORE MEMORY\n"); }
// 		printf("NBUCKETS = %d\n", NBUCKETS);
// 		printf("bitShift = %d\n", bitShift);
// 		printf("DATALENG = %d\n", DATALENG);

//     uint64_t wbuf[NBUCKETS][BUFSZ];
		uint64_t **wbuf; // mac ma maly stack
 	  wbuf = (uint64_t**)malloc(sizeof(uint64_t*) * NBUCKETS);
 		for (i = 0; i < NBUCKETS; i++) {
 			wbuf[i] = (uint64_t*)malloc(sizeof(uint64_t) * BUFSZ);
 		}
     int wbufPos[NBUCKETS];
     memset(wbufPos, 0, NBUCKETS * sizeof(int));
 
//     int wCounts[NBUCKETS]; // nmers
//     memset(wCounts, 0, NBUCKETS * sizeof(int));
//     
//     int totwCounts[NBUCKETS]; // tot nmers
//     memset(totwCounts, 0, NBUCKETS * sizeof(int));
// 
				int *dataDest = (int*)malloc(sizeof(int) * NBUCKETS);
//     int dataDest[NBUCKETS];					//pole lokacii casti dat  (part_id => destID)
//     int workerCount[NWS * NWThs];			//pocet hladin workera 
//     int cumWorkerCount[NWS * NWThs];	//zaciatocny index
//     int workerSize[NWS * NWThs];			//sum velkost dat pre workera
  
/*
    fillArray(data, DATALENG, id);
    
    MPI_Barrier(Comm_masters);
    gettimeofday(&wtStart, NULL);

///////////////////////////////////////////////////////////histogram
    gettimeofday(&preStart, NULL);
    for(i = 0; i < DATALENG; i++) {
	    wCounts[data[i] >> bitShift]++; // nmers
		}
		
//     for(i = 0; i < NBUCKETS; i++) {
// 			printf("%d: %d\n", i, wCounts[i]);
// 		}
// 		
// 		exit(2);

    /////////delenie medzi workerov
    MPI_Allreduce(&wCounts, &totwCounts, NBUCKETS, MPI_INT, MPI_SUM, Comm_masters);
    
    // partLimits ...
    uint64_t limit = (DATALENG * NMS) / (NWS * NWThs) - (( DATALENG * NMS) / NBUCKETS * 0.48);  //+0.5 z dôvodu statist. delenia priemerne 0.5 casti dalsiemu
    printf("limit = %llu\n", limit);
    int total_limit = 0, total_sum = 0;
    int sum = 0;
    int n = 0;
    int destID = 0;  //1. worker
    for(i = 0; i < NBUCKETS; i++){
        sum += totwCounts[i];
        dataDest[i] = destID;
        n++;
        
        if (sum >= limit) {
          workerCount[destID] = n;
          workerSize[destID] = sum;
          printf("destID = %d, workerCount = %d, workerSize = %d\n", destID, n, sum);      
          total_limit += limit;
          total_sum += sum;
          sum = 0;
          n = 0;
          destID = destID+1 < NWS*NWThs ? destID+1 : NWS*NWThs;          
        }        
    }
    workerCount[destID] = n;  //posledny
    workerSize[destID] = sum;
    printf("destID = %d, workerCount = %d, workerSize = %d\n", destID, n, sum);

printf("total_limit = %d, total_sum = %d\n", total_limit, total_sum);    
    cumWorkerCount[0] = 0; // partLimits[i].size ?
    for(i = 1; i < NWS * NWThs; i++){
      cumWorkerCount[i] = cumWorkerCount[i-1] + workerCount[i-1];
      printf("cumWorkerCount[%d] = %d\n", i, cumWorkerCount);
    }
  
///    MPI_Send(<#const void *buf#>, <#int count#>, <#MPI_Datatype datatype#>, <#int dest#>, <#int tag#>, <#MPI_Comm comm#>)
////////////////////////////////////////////////////////inicializacia
    if (id == 0) {
    
      for(i = 0; i < NWS; i++) {
         int initMes[NWThs*2];

         for(j = 0; j < NWThs; j++){
              initMes[j*2]   = workerCount[i*NWThs+j];      //pocet hladin (partLimits.end - partLimits.start, partLimits.count)
              initMes[j*2+1] = workerSize[i*NWThs+j];     //sumarna velkost (partLimits.size)
              printf("intMes[%d] = %d\nintMes[%d] = %d\n", j*2, workerCount[i*NWThs+j], j*2+1, workerSize[i*NWThs+j]);
         }
         printf("-----\n");
         MPI_Send((void*)initMes, 2*NWThs, MPI_INT, NMS+i, 2, MPI_COMM_WORLD);         
       }
       
       int n = 0;
       for(i = 0; i < NWS * NWThs; i++) {
           MPI_Send((void*)&(totwCounts[n]), workerCount[i], MPI_INT, i/NWThs+NMS, 3, MPI_COMM_WORLD);
           n += workerCount[i];
       }
    }
*/
    int KMER_LENGTH = 31;
       char *input_filename = "/Users/larry/Desktop/tests/raw_part_aa_withoutNs_10";
//     char *input_filename = "/Users/larry/Desktop/tests/raw_part_aa";
    char *output_path = "/Users/larry/Desktop/out/"; // TODO: nastavit na input path ak nebude zadana
    char *output_file_name = "output";
    
    int myrank = id;
    if (myrank < NMS) {
        // FILE START/END PARTITIONING
        // urcenie velkosti casti na nacitavanie (najdenie koncov riadkov)
        gettimeofday(&start, NULL);

        int parts_count = 0;
        size_t part_size;
        FileReadStartEnd *fileReadsStartsEnds;
//         printf("NMS = %d, myrank = %d, input_filename = %s\n", NMS, myrank, input_filename);
        getFileReadStartsEnds(NMS, myrank, input_filename, &fileReadsStartsEnds, &parts_count, &part_size);

        gettimeofday(&end, NULL);
        time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
        MPI_Barrier(Comm_masters);
    

    //    printf("[%d]: PARTITIONING time = %lf\n", myrank, time_spent);
    //    printf("parts_count = %d\n", parts_count);
    //    printf("part_size = %lu\n", part_size);
    //    if (myrank == 1) {
    //        for (i = 0; i < parts_count; i++) {
    //            printf("%lu: start = %lu, end = %lu\n", i, fileReadsStartsEnds[i].start, fileReadsStartsEnds[i].end);
    //        }
    //    }
        
        // KMERING (subor prechadza po castiach - kazdu cast si berie jeden proces(or)/master)
        // + NMERING
        MPI_File m_f;
        MPI_Status m_status;
        MPI_File_open(Comm_masters, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &m_f);

        
        int *nmers = (int*)malloc(sizeof(int) * nParts);
        for (i = 0; i < nParts; i++) {nmers[i] = 0;}
//         int *totNmers = (int*)malloc(sizeof(int) * nParts);
				int totNmers[nParts];
        for (i = 0; i < nParts; i++) {totNmers[i] = 0;}
       
       	int *kmersRealCounts = (int*)malloc(sizeof(int) * parts_count);
        for (i = 0; i < parts_count; i++) {kmersRealCounts[i] = -1;}
               
        int part_id;
        for (part_id = 0; part_id < parts_count; part_id++) { // niekedy musi jeden proces(or) citat aj viac casti (ak je vstupny subor prilis velky)
            int handling_rank = part_id % NMS;

            if (myrank == handling_rank) {
    //            printf("[%d]: handling_rank = %d, part_id = %d\n", myrank, handling_rank, part_id);
                printf("[%d]: %d: start = %lu, end = %lu\n", myrank, part_id, fileReadsStartsEnds[part_id].start, fileReadsStartsEnds[part_id].end);
                int bytes_to_read = (int)(fileReadsStartsEnds[part_id].end - fileReadsStartsEnds[part_id].start);
//                 printf("bytes_to_read = %d\n", bytes_to_read);
                int bytes_really_read = 0;
                char *buf = (char*)malloc(sizeof(char) * bytes_to_read);

                MPI_File_seek(m_f, fileReadsStartsEnds[part_id].start, MPI_SEEK_SET);
                MPI_File_read(m_f, buf, bytes_to_read, MPI_CHAR, &m_status);
                MPI_Get_count(&m_status, MPI_CHAR, &bytes_really_read);

                gettimeofday(&start, NULL);
                int kmers_count_approx = bytes_really_read - (bytes_really_read / KMER_LENGTH);
				        int kmers_count_real = 0;
                uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * kmers_count_approx);
                // TODO: prerobit len na n-mery, k-mery sa budu nacitavat potom ked ich bude treba
                doKmersNmers(buf, bytes_really_read, kmers, &kmers_count_real, KMER_LENGTH, nmers, nBits);
                kmersRealCounts[part_id] = kmers_count_real;
                free(kmers); // ? pouzit tie kmers niekde?
                free(buf);
    //            printf("kmers_count_approx = %d / kmers_count_real = %d\n", kmers_count_approx, kmers_count_real);
                gettimeofday(&end, NULL);
                time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
                printf("[%d]: part_id = %d, kmers = %d, K-MERING time = %lf\n", myrank, part_id, kmers_count_real, time_spent);
            }
        }
        MPI_Allreduce(nmers, &totNmers, nParts, MPI_INT, MPI_SUM, Comm_masters);
        free(nmers);

				for (i = 0; i < parts_count; i++) {
					printf("%d real kmers counts = %d\n", i, kmersRealCounts[i]);
				}

        // parts
        int parts_count_max = NWS; // pocet casti podla slave-ov
        int real_part_limits_count = 0;
        uint64_t sumNmers = 0;
        for (i = 0; i < nParts; i++) { sumNmers += totNmers[i]; }
        printf("sumNmers = %llu\n", sumNmers);
        printf("total required size for k-mers: %llu\n", sumNmers * KMER_SIZE);
        uint64_t limit = PART_SIZE_MAX;
        
        int ROUNDS = 1;
        uint64_t part_size_equal = (sumNmers * KMER_SIZE) / ( parts_count_max * ROUNDS); // aku max velkost by mala mat jedna cast
        while (part_size_equal > PART_SIZE_MAX) {
            ROUNDS++;
            part_size_equal = (sumNmers * KMER_SIZE) / (parts_count_max * ROUNDS);
        }
//         real_part_limits_count = ROUNDS * parts_count;

        if (myrank == 0) {
            printf("MAX   PART SIZE = %d\n", PART_SIZE_MAX);
            printf("EQUAL PART SIZE = %llu\n", part_size_equal);
            printf("ROUNDS = %d\n", ROUNDS); // uz nepouzivame
//             printf("REAL PARTS COUNT (ROUNDS * NWThs) = %d\n", real_part_limits_count);
        }
        
        // jeden slave by mal zvladnut spravit jeden z partLimits
        PartLimit *partLimits = (PartLimit*)malloc(sizeof(PartLimit) * nParts); // viac casti nemozeme mat
        int indexToPartId[nParts];
        
        int curr_part = 0;          // id casti
        int num_parts_current = 0;  // kolko casti obsahuje (n)
        int part_size_curr = 0;     // aku ma tato cast velkost
        int part_count_curr = 0;    // kolko casti obsahuje
        partLimits[curr_part].id = curr_part;
        partLimits[curr_part].start = 0;
        partLimits[curr_part].end = 0;
        indexToPartId[0] = partLimits[curr_part].id;
        for (i = 0; i < nParts; i++) {
            part_count_curr += totNmers[i];
            part_size_curr += totNmers[i] * KMER_SIZE;
            num_parts_current++;
		        indexToPartId[i] = partLimits[curr_part].id;

            if (part_size_curr > part_size_equal || i == nParts - 1) {
                partLimits[curr_part].end = (int)i;
                partLimits[curr_part].size = part_size_curr;
                partLimits[curr_part].count = part_count_curr;
                partLimits[curr_part].num_parts = num_parts_current;
                
                curr_part++;
                partLimits[curr_part].id = curr_part; // == destID ?
                partLimits[curr_part].start = (int)i + 1;
                partLimits[curr_part].end = (int)i + 1;
                part_size_curr = 0;
                num_parts_current = 0;
                part_count_curr = 0;
				        indexToPartId[i] = partLimits[curr_part].id;
                
                real_part_limits_count++;
            }
        }
				printf("REAL PARTS COUNT = %d\n", real_part_limits_count);

        if (myrank == 0) {
            for (i = 0; i < real_part_limits_count; i++) {
                printf("part_id = %d, start = %d, end = %d, num_parts = %d, count = %llu, size = %llu\n", partLimits[i].id, partLimits[i].start, partLimits[i].end, partLimits[i].num_parts, partLimits[i].count, partLimits[i].size);
            }
        }
        MPI_Barrier(Comm_masters);
        
        if (myrank == 0) {
            int round;
            
            for (round = 0; round < ROUNDS; round++) {
              int initMes[NWThs * 2];
              int initMes2[NWS][NWThs * 2];
							int dataDestWorker[nParts];
							int dataDestThread[nParts];
							for (i = 0; i < nParts; i++) {
								dataDestWorker[i] = -1;
								dataDestThread[i] = -1;
							}
            	
                for (i = 0; i < NWS; i++) {
                    int previousJ = 0;
                    int roundNWSId = (round * ROUNDS) + i;

										printf("[%d]: ROUND = %d\n", myrank, round);
										printf("FROM %d TO %d (num_parts = %d)\n", partLimits[roundNWSId].start, partLimits[roundNWSId].end, partLimits[roundNWSId].num_parts);

										int slave_num_parts = -1;
										int slave_counts = -1;

										// kazdy part sa rozdeli podla poctu slave-ov
// 										if (NWS > 1) {
// 										
// 											if (i == 0) { // priratat aj zvysok (jedenkrat)
// 													slave_num_parts   = partLimits[roundNWSId].num_parts % NWS == 0 ? partLimits[roundNWSId].num_parts / NWS : partLimits[roundNWSId].num_parts / NWS + partLimits[roundNWSId].num_parts % NWS;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
// 													slave_counts = partLimits[roundNWSId].count % NWS == 0 ? partLimits[roundNWSId].count / NWS : partLimits[roundNWSId].count / NWS + partLimits[roundNWSId].count % NWS;     //sumarna velkost (partLimits.size)
// 											} else {
// 													slave_num_parts   = partLimits[roundNWSId].num_parts / 2;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
// 													slave_counts = partLimits[roundNWSId].count / NWS;     //sumarna velkost (partLimits.size)
// 											}
// 										} else {
											slave_num_parts = partLimits[roundNWSId].num_parts;
											slave_counts = partLimits[roundNWSId].count;
// 										}


                    for (j = 0; j < NWThs; j++) {
// 	                      printf("worker = %d, thread = %d\n", i, j);
                        // kazdy part sa rozdeli posla poctu threadov
                        
                        if (NWThs > 1) {

													if (j == 0) { // priratat aj zvysok (jedenkrat)
															initMes[j*2]   = slave_num_parts % NWThs == 0 ? slave_num_parts / NWThs : slave_num_parts / NWThs + slave_num_parts % NWThs;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
															initMes[j*2+1] = slave_counts % NWThs == 0 ? slave_counts / NWThs : slave_counts / NWThs + slave_counts % NWThs;     //sumarna velkost (partLimits.size)
													} else {
															initMes[j*2]   = slave_num_parts / NWThs;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
															initMes[j*2+1] = slave_counts / NWThs;     //sumarna velkost (partLimits.size)
													}
// 													if (j == 0) { // priratat aj zvysok (jedenkrat)
// 															initMes[j*2]   = partLimits[roundNWSId].num_parts % NWThs == 0 ? partLimits[roundNWSId].num_parts / NWThs : partLimits[roundNWSId].num_parts / NWThs + partLimits[roundNWSId].num_parts % NWThs;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
// 															initMes[j*2+1] = partLimits[roundNWSId].count % NWThs == 0 ? partLimits[roundNWSId].count / NWThs : partLimits[roundNWSId].count / NWThs + partLimits[roundNWSId].count % NWThs;     //sumarna velkost (partLimits.size)
// 													} else {
// 															initMes[j*2]   = partLimits[roundNWSId].num_parts / NWThs;      //pocet hladin (partLimits.end - partLimits.start, partLimits.num_parts)
// 															initMes[j*2+1] = partLimits[roundNWSId].count / NWThs;     //sumarna velkost (partLimits.size)
// 													}
                        } else {
                        	printf("\n NWThs < 1 !!!!!!\n\n");
   												initMes[j*2]   = slave_num_parts;
		  										initMes[j*2+1] = slave_counts;
                        }
                        initMes2[i][j*2] = initMes[j*2];
                        initMes2[i][j*2+1] = initMes[j*2+1];
                        
                        int q;
                        for (q = partLimits[roundNWSId].start + previousJ; q < partLimits[roundNWSId].start + previousJ + initMes[j*2]; q++) { // ktory thread sa bude o toto starat
                        	dataDestThread[q] = j;
                        }
                        previousJ = initMes[j*2];
                        
                        printf("intMes[%d] = %d\nintMes[%d] = %d\n", j*2, initMes[j*2], j*2+1, initMes[j*2+1]);
                    }
                    int q;
                    for (q = partLimits[roundNWSId].start; q < partLimits[roundNWSId].start + partLimits[roundNWSId].num_parts; q++) { // ktory worker sa bude o toto starat
                    	dataDestWorker[q] = i;
                    }
//                     printf("---next NWS--\n");
//                     printf("Master %d sending TAG 2 to %d\n", myrank, NMS+(int)i);
//                     printf("%d %d %d %d\n", initMes[0], initMes[1], initMes[2], initMes[3]);
                    MPI_Send((void*)&initMes, 2 * NWThs, MPI_INT, NMS + (int)i, 2, MPI_COMM_WORLD);
                }

// 								for (i = 0; i < nParts; i++) {
// 									printf("%3d: worker: %d | thread: %d\n", i, dataDestWorker[i], dataDestThread[i]);
// 								}
// 								MPI_Finalize();
// 								exit(0);

//                int n = 0;
//                for (i = 0; i < NWS * NWThs; i++) { // posle sa pocet nParts a ich celkova pocetnost (totNmers)
// //                    printf("Master %d sending TAG 3 to %d\n", myrank, (int)i / NWThs + NMS);
// //                    printf("totNmers[%d] = %d, initMes[%d] = %d\n", n, totNmers[n], i, initMes[i*NWThs]);
//                    MPI_Send((void*)(&totNmers[n]), initMes[i*NWThs], MPI_INT, (int)i / NWThs + NMS, 3, MPI_COMM_WORLD);
//                    printf("n = %d\n", n);
//                    n += partLimits[round].num_parts - 1;
//                }
									int n = 0;
									for (i = 0; i < NWS; i++) {

										for (j = 0; j < NWThs; j++) {
 							          MPI_Send((void*)&(totNmers[n]), initMes2[i][j*2], MPI_INT, NMS + (i*NWS+j)/NWThs, 3, MPI_COMM_WORLD);
                        n += initMes2[i][j*2];
										}
									}

//                int *hlp;
//                MPI_Status recvStatus;
//                MPI_Recv(&hlp, 1, MPI_INT, MPI_ANY_SOURCE, 444, MPI_COMM_WORLD, &recvStatus); 
//                printf("Master: 444\n");
							
// 							int roundStartEnd[2]; // poslat vsetkym masterom start/end ktore maju nacitat zo suboru
// 							roundStartEnd[0] = partLimits[round].start;
// 							roundStartEnd[1] = partLimits[round].end;

							for (i = 1; i < NMS; i++) {
							  // poslat masterom info o tom ktore k-tice kam pojdu (ktoremu workerovi a ktoremu vlaknu)
							  MPI_Send((void*)&dataDestWorker, nParts, MPI_INT, i, 111, Comm_masters);
							  MPI_Send((void*)&dataDestThread, nParts, MPI_INT, i, 222, Comm_masters);
							  // aby zacali aj ostatni masteri prechadzat subor a vyhladavat dane k-tice
								MPI_Send((void*)&round, 1, MPI_INT, i, 555, Comm_masters);
							}
							
							// rank 0 bude tiez robit svoju cast (spolu s ostatnymi mastermi)
							// nacitat k-tice
							int part_limit_id;
 							for (part_limit_id = 0; part_limit_id < real_part_limits_count; part_limit_id++) {
								printf("vyberam k-tice od %d do %d (celkovy pocet = %d)\n", partLimits[part_limit_id].start, partLimits[part_limit_id].end, partLimits[part_limit_id].count);
								for (part_id = 0; part_id < parts_count; part_id++) { // ked je velky subor, tak jeden proces(or) musi robit viac casti...
									int handling_rank = part_id % NMS;
								
									if (handling_rank == myrank) {
										uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * partLimits[part_limit_id].count); // partLimits[part_limit_id].count == pocet k-tic
										int kmers_count_real = 0;
										// void getKmers(MPI_File m_f, FileReadStartEnd *fileReadsStartsEnds, int part_id, int start_level, int end_level, uint64_t *kmers, int KMER_LENGTH, int nBits) {
										getKmers(m_f, fileReadsStartsEnds, part_id, partLimits[part_limit_id].start, partLimits[part_limit_id].end, kmers, KMER_LENGTH, nBits, &kmers_count_real);

// 										for (i = 0; i < kmers_count_real; i++) {
// 											printf("%llu\n", kmers[i]);
// 										}

										// ROZPOSIELANIE WORKEROM

										int bitShift = 64 - nBits;
										for (i = 0; i < kmers_count_real; i++) {
											int bInx = kmers[i] >> bitShift;
											wbuf[bInx][wbufPos[bInx]++] = kmers[i];
										
											if (wbufPos[bInx] == BUFSZ) {
													int dataLMes[3];
													int worker = dataDestWorker[bInx];
													int thread = dataDestThread[bInx];
													int inx = partLimits[indexToPartId[bInx]].start;
													dataLMes[0] = BUFSZ;            //size
													dataLMes[1] = thread;       //dest thread
													dataLMes[2] = inx;          //dest inx
													MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
													MPI_Send((void*)wbuf[bInx], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
													wbufPos[bInx] = 0;
											}
										}
									
										// odoslanie zvyskov ktore nedovrsili BUFSZ
										for (i = 0; i < NBUCKETS; i++) { // NBUCKETS ~ nParts
		
												if (wbufPos[i] == 0) {   //nezasli falosnu end mes
													continue;
												}
												int dataLMes[3];
												int worker = dataDestWorker[i];
												int thread = dataDestThread[i];
												int inx = partLimits[indexToPartId[i]].start;
												dataLMes[0] = wbufPos[i];            //size
												dataLMes[1] =  thread;       //dest thread
												dataLMes[2] = inx;          //dest inx
												MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
												MPI_Send((void*)wbuf[i], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
										}
		
										// start?							
										for (i = 0; i < NWS; i++) {
											for (j = 0; j < NWThs; j++) {
												 int dataLMes[3];
												 dataLMes[0] = 0;
												 dataLMes[1] = j;
												 MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+i, 0, MPI_COMM_WORLD);   //end message
											}
										}
									
									
									
										free(kmers);
									}
								}
      		  	}



            }

        } else { // end if myrank == 0 // ostatni mastri
					
					while (1) {
						MPI_Status recvStatus;
						int round;
						int dataDestWorker[nParts];
						int dataDestThread[nParts];
						MPI_Recv(&dataDestWorker, nParts, MPI_INT, MPI_ANY_SOURCE, 111, Comm_masters, &recvStatus);
						MPI_Recv(&dataDestThread, nParts, MPI_INT, MPI_ANY_SOURCE, 222, Comm_masters, &recvStatus);
						MPI_Recv(&round, 1, MPI_INT, MPI_ANY_SOURCE, 555, Comm_masters, &recvStatus);
						printf("Master received to do round: %d\n", round);

							// nacitat k-tice
							int part_limit_id;
 							for (part_limit_id = 0; part_limit_id < real_part_limits_count; part_limit_id++) {
								printf("vyberam k-tice od %d do %d (celkovy pocet = %d)\n", partLimits[part_limit_id].start, partLimits[part_limit_id].end, partLimits[part_limit_id].count);
								for (part_id = 0; part_id < parts_count; part_id++) { // ked je velky subor, tak jeden proces(or) musi robit viac casti...
									int handling_rank = part_id % NMS;
								
									if (handling_rank == myrank) {
										uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * partLimits[part_limit_id].count); // partLimits[part_limit_id].count == pocet k-tic
										int kmers_count_real = 0;
										// void getKmers(MPI_File m_f, FileReadStartEnd *fileReadsStartsEnds, int part_id, int start_level, int end_level, uint64_t *kmers, int KMER_LENGTH, int nBits) {
										getKmers(m_f, fileReadsStartsEnds, part_id, partLimits[part_limit_id].start, partLimits[part_limit_id].end, kmers, KMER_LENGTH, nBits, &kmers_count_real);

// 										for (i = 0; i < kmers_count_real; i++) {
// 											printf("%llu\n", kmers[i]);
// 										}

										// ROZPOSIELANIE WORKEROM

										int bitShift = 64 - nBits;
										for (i = 0; i < kmers_count_real; i++) {
											int bInx = kmers[i] >> bitShift;
											wbuf[bInx][wbufPos[bInx]++] = kmers[i];
										
											if (wbufPos[bInx] == BUFSZ) {
													int dataLMes[3];
													int worker = dataDestWorker[bInx];
													int thread = dataDestThread[bInx];
													int inx = partLimits[indexToPartId[bInx]].start;
													dataLMes[0] = BUFSZ;            //size
													dataLMes[1] = thread;       //dest thread
													dataLMes[2] = inx;          //dest inx
													MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
													MPI_Send((void*)wbuf[bInx], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
													wbufPos[bInx] = 0;
											}
										}
									
										// odoslanie zvyskov ktore nedovrsili BUFSZ
										for (i = 0; i < NBUCKETS; i++) { // NBUCKETS ~ nParts
		
												if (wbufPos[i] == 0) {   //nezasli falosnu end mes
													continue;
												}
												int dataLMes[3];
												int worker = dataDestWorker[i];
												int thread = dataDestThread[i];
												int inx = partLimits[indexToPartId[i]].start;
												dataLMes[0] = wbufPos[i];            //size
												dataLMes[1] =  thread;       //dest thread
												dataLMes[2] = inx;          //dest inx
												MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
												MPI_Send((void*)wbuf[i], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
										}
		
										// start?							
										for (i = 0; i < NWS; i++) {
											for (j = 0; j < NWThs; j++) {
												 int dataLMes[3];
												 dataLMes[0] = 0;
												 dataLMes[1] = j;
												 MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+i, 0, MPI_COMM_WORLD);   //end message
											}
										}
									
										free(kmers);
									}
								}
      		  	}
						
					}  // /while 1
      	}
        
        
        
        free(partLimits);
        MPI_File_close(&m_f);
		} else { // id > NMS
			// TODO:
			while(1) {
				// spustit novy worker code pre kazde "nove kolo"
				workerCode(id); // odkomentoavt a spravit moznost opakovania volania worker kodu
				
				int new_round;
				MPI_Status recvStatus;
				MPI_Recv(&new_round, 1, MPI_INT, MPI_ANY_SOURCE, 777, MPI_COMM_WORLD, &recvStatus);
			}
		}

// exit(1);
/*
/////////////////////////////////////////////////////////odosielanie
    gettimeofday(&start, NULL);

    for(i = 0; i < DATALENG; i++) {         //seriove, inak kontraproduktivne
        int bInx = data[i] >> bitShift;

        wbuf[bInx][wbufPos[bInx]++] = data[i];
        
        if(wbufPos[bInx] == BUFSZ) {
            int dataLMes[3];
            int worker = dataDest[bInx] / NWThs;
            int thread = dataDest[bInx] % NWThs;
            int inx = bInx - cumWorkerCount[dataDest[bInx]];
            dataLMes[0] = BUFSZ;            //size
            dataLMes[1] =  thread;       //dest thread
            dataLMes[2] = inx;          //dest inx
            MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
            MPI_Send((void*)wbuf[bInx], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
            wbufPos[bInx] = 0;
        }
    }
 
    for(i = 0; i < NBUCKETS; i++) {
    
        if (wbufPos[i] == 0) {   //nezasli falosnu end mes
          continue;
        }
        int dataLMes[3];
        int worker = dataDest[i] / NWThs;
        int thread = dataDest[i] % NWThs;
        int inx = i - cumWorkerCount[dataDest[i]];
        dataLMes[0] = wbufPos[i];            //size
        dataLMes[1] =  thread;       //dest thread
        dataLMes[2] = inx;          //dest inx
        MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+worker, 0, MPI_COMM_WORLD);
        MPI_Send((void*)wbuf[i], dataLMes[0], MPI_UNSIGNED_LONG_LONG, NMS+worker, 1, MPI_COMM_WORLD);
    }
    
    for(i=0;i<NWS;i++){
      for(j=0;j<NWThs;j++){
         int dataLMes[3];
         dataLMes[0] = 0;
         dataLMes[1] = j;
         MPI_Send((void*)dataLMes, 3, MPI_LONG, NMS+i, 0, MPI_COMM_WORLD);   //end message
      }
    }

    gettimeofday(&end, NULL);
    timeHist = (start.tv_sec - preStart.tv_sec) * 1000 + (start.tv_usec - preStart.tv_usec)/1000.0;
    timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
    
    free(data);
    int length;
    char hostID[100];
    MPI_Get_processor_name(hostID, &length);  
    
    printf("%d, %s: preproc %lf, odosielanie %lf\n", id, hostID, timeHist, timeSpent);
    
//     MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&wtEnd, NULL);
    if(id == 0){
         timeSpent = (wtEnd.tv_sec - wtStart.tv_sec) * 1000 + (wtEnd.tv_usec - wtStart.tv_usec)/1000.0;
         printf("WALLCLOCKTIME %lf\n", timeSpent);
    }
    */
}


int main(int argc,char *argv[]) {

    BASES[65]=0; BASES[67]=1; BASES[71]=2; BASES[84]=3;
    BASES_2[65]=0<<2; BASES_2[67]=1<<2; BASES_2[71]=2<<2; BASES_2[84]=3<<2;
    
		
		DATALENG = 10000000;//599890416 / 8; // aa 6 //599890416

    int myid, numprocs;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

		NWS = (int)(numprocs / 2);            //pocet workerov
		NWThs = 2;          //pocet vlakien workera
		NMS = (int)(numprocs / 2);            //pocet masterov
		NMThs = 1;          //pocet vlakien mastera //OPENMPI - BUSY WAITING, ine ako 1 nic neriesi
   
    struct timeval ts;
    gettimeofday(&ts, NULL);
    srand(ts.tv_usec * (ts.tv_sec+myid));
    
    MPI_Group orig_group, new_group;  
    MPI_Comm dup_comm_world;        
    MPI_Comm_dup( MPI_COMM_WORLD, &dup_comm_world );               
    MPI_Comm_group(dup_comm_world, &orig_group);         
    int mRanks[NMS]; 
    int i;
    
    for (i = 0; i < NMS; i++) {
    	mRanks[i] = i;
		}
    MPI_Group_incl(orig_group, NMS, mRanks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &Comm_masters);

		MPI_Barrier(MPI_COMM_WORLD);
    if (myid >= NMS+NWS){
        MPI_Finalize();
        exit(0);
    }
    
//     if (myid < NMS) {
// 				printf("master id = %d\n", myid);
        masterCode(myid);
//     } else {
//         printf("worker id = %d\n", myid);
//         workerCode(myid);
//    }

    MPI_Finalize();
    exit(0);
}
