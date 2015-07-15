#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// kmers extraction
#define MAX_MEMORY_PROCESSOR_FILE_LOAD 1073741824
#define RESERVE 1024
#define PART_SIZE_MAX 1073741824
#define KMER_SIZE 8 // 8B

// sorting
#define maxLgNBs 12

char BASES[255];
char BASES_2[255];

typedef struct PartLimit {
    int start;
    int end;
    int num_parts;
    uint64_t size;
    size_t count;
} PartLimit;

typedef struct KeyVal {
    uint64_t kmer;
    int count;
} KeyVal;

int myPow(int n, int nn) {
    int i, tmp = n;
    for (i = 1; i < nn; i++) {
        tmp = tmp * n;
    }
    return nn == 0 ? 1 : tmp;
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
    int i, j;
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


void bucketSortSerial(uint64_t *data, uint64_t *aux, uint64_t dataCount, int _baseBitShift, int lgNBUCKETS, short beNested) { //vstup v aux, vystup v data   // dataCount - pocet prvkov, _bbs - pocet zhora rovnakych bitov, lgNB - log. pocetu vedierok
    int i;
    
    int NBUCKETS = (1 << lgNBUCKETS);
    int bitShift = 64 - lgNBUCKETS;
    
    int bCounts[NBUCKETS];       memset(bCounts, 0, NBUCKETS * sizeof(int));
    int cumBCounts[NBUCKETS];    memset(cumBCounts, 0, NBUCKETS * sizeof(int));//aka bucket start
    
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
        if (beNested == 0) {
            for(i=0; i<NBUCKETS; i++) {
                quick_sort(&data[cumBCounts[i]], (int*)&aux[cumBCounts[i]], bCounts[i]);
            }
        } else {
            for(i=0; i<NBUCKETS; i++) {
                int n = floor(log((double)bCounts[i] / 25) / log(2) );
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


// koniec sekvencie - bud \n (pre raw) alebo \n'>|@' pre fasta
size_t getEndOfSequence(int INPUT_FILE_TYPE, char *buf, size_t size_to_load, size_t reserve) { // 0: raw, 1: fasta
    size_t i;
    size_t end_offset = 0;
    for (i = size_to_load - reserve; i < size_to_load; i++) {
        end_offset++;
        if (buf[i] == '\n') {
            return end_offset;
        }
    }
    return 0;
}

size_t getEndOfSequenceFromForward(int INPUT_FILE_TYPE, char *buf, size_t from, size_t reserve) { // 0: raw, 1: fasta
    size_t i;
    size_t end_offset = 0;
    for (i = from; i < from + reserve; i++) {
        end_offset++;
        if (buf[i] == '\n') {
            return end_offset;
        }
    }
    return 0;
}

// funkcia zisti kde su konce riadkov ap odla toho nastavi jednotlive hodnoty pre procesy - ktory proces ma odkial-pokial citat suboru
void getPartsStartEnd(int myrank, int npes, char *filename, size_t *process_part_read_start, size_t *process_part_read_end) {
    
    FILE *f = fopen(filename, "r");
    
    fseek(f, 0, SEEK_END);
    size_t f_size = ftell(f);
    fseek(f, 0, SEEK_SET);
//    printf("filesize = %lu\n", f_size);
    
    size_t process_part_size = (f_size % npes == 0) ? f_size / npes : f_size / npes + f_size % npes; // aku velku cast z celkovej velkosti suboru bude tento proces spracovavat
    //    printf("process_part_size = %lu\n", process_part_size);
    
    *process_part_read_start = myrank * process_part_size;
    *process_part_read_end = (myrank+1) * process_part_size;
    
    size_t process_move_start_by[npes]; // o kolko posunut start pri citani zo suboru
    size_t process_move_end_by[npes]; // o kolko posunut end pri citani zo suboru
    size_t move_start_by[npes];
    size_t move_end_by[npes];
    
    int i;
    for (i = 0; i < npes; i++) {
        process_move_start_by[i] = 0;
        process_move_end_by[i] = 0;
        move_start_by[i] = 0;
        move_end_by[i] = 0;
    }
    //    printf("[%d]: process_part_read_start = %lu => process_part_read_end = %lu\n", myrank, process_part_read_start, process_part_read_end);
    // skocit na koniec a zobrat nejaku rezervu dopredu
    size_t reserve = RESERVE;
    size_t size_to_load = sizeof(char) * reserve;
    char *buf = (char*)malloc(size_to_load);
    
    if (myrank == 0) {
        // najst novy end
        fseek(f, *process_part_read_end, SEEK_SET);
        size_t loaded_size = fread(buf, sizeof(char), size_to_load, f);
        
        size_t end_offset = getEndOfSequence(1, buf, size_to_load, reserve); // najst najblizsi koniec riadka (smerom dopredu)
        process_move_end_by[myrank] = end_offset;
        memset(buf, '\0', size_to_load);
    } else { // myrank != 0
        // najst novy start
        fseek(f, *process_part_read_start, SEEK_SET);
        size_t loaded_size = fread(buf, sizeof(char), size_to_load, f);
        size_t start_offset = getEndOfSequence(1, buf, size_to_load, reserve); // najst najblizsi koniec riadka (smerom dopredu)
        process_move_start_by[myrank] = start_offset;
        memset(buf, '\0', size_to_load);
        
        // najst novy end
        fseek(f, *process_part_read_end, SEEK_SET);
        loaded_size = fread(buf, sizeof(char), size_to_load, f);
        size_t end_offset = getEndOfSequence(1, buf, size_to_load, reserve); // najst najblizsi koniec riadka (smerom dopredu)
        process_move_end_by[myrank] = end_offset;
        memset(buf, '\0', size_to_load);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&process_move_start_by, &move_start_by, npes, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&process_move_end_by, &move_end_by, npes, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
    
    // fix starts and ends
    *process_part_read_end += move_end_by[myrank];
    if (myrank  > 0) {
        
        if (*process_part_read_end > f_size) { // aby posledny proces nepreskocil velkost suboru
            *process_part_read_end = f_size;
        }
        *process_part_read_start += move_end_by[myrank-1]; // posunut start o end predosleho
    }
    free(buf);
    fclose(f);
}

// vytiahnne z bufferu raw readov n-tice (a pripocita ich do nmers)
void getNmersFromRaw(char *buf, size_t buf_size, int KMER_LENGTH, int nBases, int *nmers) {
    int nBits = nBases * 2; // => 2^nBits roznych casti, podla vrchnych bitov
    int nmask = 0;
    size_t i;
    for (i = 0; i < nBits; i++) { // nastavenie "nulovacej" masky -> 00000....1111...1
        nmask |= (1 << i);
    }
    
    int nmer_length_counter = 0;
    int nmer = 0;
    int nmers_counter = 0;
    for (i = 0; i < buf_size; i++) {
        
        if (buf[i] != '\n' && buf[i] != '\0' && buf[i] != EOF && buf[i] != 'N') {
            
            if (nmer_length_counter < nBases) {
                nmer |= (BASES[buf[i]] << (nBits - nmer_length_counter * 2 - 2));
                nmer_length_counter++;
                if (nmer_length_counter == nBases) {
                    nmers[nmer]++;
                    nmers_counter++;
                }
            } else {
                nmer = ((nmer << 2) & nmask) + BASES[buf[i]];
                nmers[nmer]++;
                nmers_counter++;
            }
        } else {
            
            if (nmers_counter > 0) { // musim odstranit poslednych |K|-|N| n-tic
                size_t j;
                for (j = i-1; j > i-KMER_LENGTH+nBases-1; j--) {
                    nmers[nmer]--;
                    nmer = (nmer >> 2) + (BASES[buf[j-nBases]] << (nBits - 2));
                }
            }
            nmer = 0;
            nmer_length_counter = 0;
            nmers_counter = 0;
        }
    }
}

void getPartNmers(int myrank, MPI_File m_f, size_t process_part_read_start, size_t process_part_read_end, int KMER_LENGTH, int nBases, int *nmers) {
    size_t process_part_read_current = process_part_read_start;
    //    FILE *f = fopen(input_filename, "r");
    MPI_Status m_status;
    MPI_File_seek(m_f, process_part_read_current, MPI_SEEK_SET);
//    fseek(f, process_part_read_current, SEEK_SET);
    char *buf = (char*)malloc(sizeof(char) * MAX_MEMORY_PROCESSOR_FILE_LOAD + RESERVE);
    
    size_t reserve_size = 0;
    size_t reserve_size_old = reserve_size;
    while (process_part_read_current < process_part_read_end) {
        size_t read_next_bytes;
        if (process_part_read_current + MAX_MEMORY_PROCESSOR_FILE_LOAD > process_part_read_end) { // neprekrocit koniec max. pamatou
            read_next_bytes = process_part_read_end - process_part_read_current;
        } else {
            read_next_bytes = MAX_MEMORY_PROCESSOR_FILE_LOAD;
        }
//        size_t loaded_size = fread(buf, sizeof(char), read_next_bytes + RESERVE, f);
        read_next_bytes = read_next_bytes + RESERVE;
        if (process_part_read_current + read_next_bytes > process_part_read_end) { // nepreskocit koniec rezervou
            read_next_bytes = process_part_read_end - process_part_read_current;
        }
        size_t loaded_size = 0;
        MPI_File_read(m_f, buf, read_next_bytes, MPI_CHAR, &m_status);
        MPI_Get_count(&m_status, MPI_CHAR, &loaded_size);
        reserve_size_old = reserve_size;
        reserve_size = getEndOfSequenceFromForward(1, buf, read_next_bytes, RESERVE); // koniec sekvencie od (from) smerom dopredu (forward)
//        fseek(f, process_part_read_current + read_next_bytes + reserve_size, SEEK_SET);
        MPI_File_seek(m_f, process_part_read_current + read_next_bytes + reserve_size, MPI_SEEK_SET);
        memset(&buf[read_next_bytes + reserve_size], '\0', 1);
//        printf("[%d]: %lu -> %lu: %s", myrank, process_part_read_current, process_part_read_current + read_next_bytes + reserve_size, buf);
//        printf("%s", buf);
        getNmersFromRaw(buf, read_next_bytes + reserve_size, KMER_LENGTH, nBases, nmers);
        memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + RESERVE);

        if (process_part_read_current + read_next_bytes + reserve_size >= process_part_read_end) {
            break;
        }
        process_part_read_current += read_next_bytes + reserve_size;
    }
//    fclose(f);
    free(buf);
}

void getKmersFromRaw(char *buf, size_t buf_size, int KMER_LENGTH, int nBases, int start_level, int end_level, uint64_t *kmers, int *kmers_count) {
    int nBits = nBases * 2; // => 2^nBits roznych casti, podla vrchnych bitov
    int nmask = 0;
    size_t i;
    for (i = 0; i < nBits; i++) { // nastavenie "nulovacej" masky -> 00000....1111...1
        nmask |= (1 << i);
    }
    
    uint64_t kmer = 0;
    uint64_t nmer_value = 0;
    int kmer_length_counter = 0;
//    printf("%s", buf);
    for (i = 0; i < buf_size; i++) {
        
        if (buf[i] != '\n' && buf[i] != '\0' && buf[i] != EOF && buf[i] != 'N') {
            
            if (kmer_length_counter < KMER_LENGTH) {
                kmer |= (((uint64_t)BASES[buf[i]]) << (62 - kmer_length_counter * 2));
                kmer_length_counter++;
                if (kmer_length_counter == KMER_LENGTH) {
                    nmer_value = (kmer >> (64 - nBits)); // ziskat hodnotu prvych nBits tejto k-tice
                    
                    if (nmer_value >= start_level && nmer_value <= end_level) {
                        kmers[(*kmers_count)] = kmer;
                        (*kmers_count)++;
//                        printf("first kmer = %llu\n", kmer);
                    }
                }
            } else {
                kmer = (uint64_t)((kmer << 2) | ((uint64_t)BASES_2[buf[i]])); // << 2
                nmer_value = (kmer >> (64 - nBits)); // ziskat hodnotu prvych nBits tejto k-tice
                
                if (nmer_value >= start_level && nmer_value <= end_level) {
                    kmers[(*kmers_count)] = kmer;
                    (*kmers_count)++;
//                    printf("kmer = %llu\n", kmer);
                }
            }
        } else {
            kmer = 0;
            kmer_length_counter = 0;
        }
    }
//    printf("kmers_count = %d\n", *kmers_count);
}

void getPartKmers(int myrank, MPI_File m_f, int start_level, int end_level, int KMER_LENGTH, int nBases, uint64_t *kmers, int *kmers_count, size_t *parts_read_start, size_t *parts_read_end) {
    size_t process_part_read_current = 0;
    MPI_Status m_status;
    MPI_Offset m_filesize;
    MPI_File_get_size(m_f, &m_filesize);
    MPI_File_seek(m_f, 0, MPI_SEEK_SET);
    char *buf = (char*)malloc(sizeof(char) * MAX_MEMORY_PROCESSOR_FILE_LOAD + RESERVE);
    size_t reserve_size = 0;
    size_t reserve_size_old = reserve_size;
    while (process_part_read_current < m_filesize) {
        size_t read_next_bytes;
        if (process_part_read_current + MAX_MEMORY_PROCESSOR_FILE_LOAD > m_filesize) { // neprekrocit koniec max. pamatou
            read_next_bytes = m_filesize - process_part_read_current;
        } else {
            read_next_bytes = MAX_MEMORY_PROCESSOR_FILE_LOAD;
        }
        //        size_t loaded_size = fread(buf, sizeof(char), read_next_bytes + RESERVE, f);
        read_next_bytes = read_next_bytes + RESERVE;
        if (process_part_read_current + read_next_bytes > m_filesize) { // nepreskocit koniec rezervou
            read_next_bytes = m_filesize - process_part_read_current;
        }
        size_t loaded_size = 0;
        MPI_File_read(m_f, buf, read_next_bytes, MPI_CHAR, &m_status);
        MPI_Get_count(&m_status, MPI_CHAR, &loaded_size);
//        printf("%s", buf);
        reserve_size_old = reserve_size;
        reserve_size = getEndOfSequenceFromForward(1, buf, read_next_bytes, RESERVE); // koniec sekvencie od (from) smerom dopredu (forward)
        MPI_File_seek(m_f, process_part_read_current + read_next_bytes + reserve_size, MPI_SEEK_SET);
        memset(&buf[read_next_bytes + reserve_size], '\0', 1);
        getKmersFromRaw(buf, read_next_bytes + reserve_size, KMER_LENGTH, nBases, start_level, end_level, kmers, kmers_count);
        memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + RESERVE);
        
        if (process_part_read_current + read_next_bytes + reserve_size >= m_filesize) {
            break;
        }
        process_part_read_current += read_next_bytes + reserve_size;
    }
    free(buf);
}

void countPartsLimits(int myrank, int npes, int nParts, int *totNmers, PartLimit *partLimits, int *num_parts) {
    // bolo by dobre aby casti boli rovnako velke, casti je 'parts_count_max'
    int parts_count_max = npes; // podla poctu procesorov (uzlov, ...)
    uint64_t kmers_count_total = 0;

    size_t i;
    for (i = 0; i < nParts; i++) {
        kmers_count_total += totNmers[i];
    }

    int rounds = 1;
    uint64_t part_size_equal = (kmers_count_total * KMER_SIZE) / (parts_count_max * rounds); // aku max velkost by mala mat jedna cast
    while (part_size_equal > PART_SIZE_MAX) {
        rounds++;
        part_size_equal = (kmers_count_total * KMER_SIZE) / (parts_count_max * rounds);
    }
    if (myrank == 0) {
        printf("MAX   PART SIZE = %d\n", PART_SIZE_MAX);
        printf("EQUAL PART SIZE = %llu\n\n", part_size_equal);
        //        printf("ROUNDS = %d\n", rounds); // uz nepouzivame
    }
    // rozdelit na rovnomerne casti
//    int num_parts = 0;
    int num_parts_current = 0;
    uint64_t part_size_curr = 0;
    uint64_t part_count_curr = 0;
    partLimits[*num_parts].start = 0;
    partLimits[*num_parts].end = 0;
    for (i = 0; i < nParts; i++) {
        part_count_curr += totNmers[i];
        part_size_curr += totNmers[i] * KMER_SIZE;
        num_parts_current++;
        
        if (part_size_curr > part_size_equal || i == nParts - 1) {
            partLimits[*num_parts].end = (int)i;
            partLimits[*num_parts].size = part_size_curr;
            partLimits[*num_parts].count = part_count_curr;
            partLimits[*num_parts].num_parts = num_parts_current;
            //            printf("[%d] set id = %d, start = %d, end = %d\n", myrank, num_parts, partLimits[num_parts].start, partLimits[num_parts].end);
            (*num_parts)++;
            partLimits[*num_parts].start = (int)i + 1;
            partLimits[*num_parts].end = (int)i + 1;
            part_size_curr = 0;
            num_parts_current = 0;
            part_count_curr = 0;
        }
    }
}

int main(int argc, char * argv[]) {
    struct timeval start, end, start_total, end_total;
    gettimeofday(&start_total, NULL);
    double time_spent;
    
    BASES[65]=0; BASES[67]=1; BASES[71]=2; BASES[84]=3;
    BASES_2[65]=0<<2; BASES_2[67]=1<<2; BASES_2[71]=2<<2; BASES_2[84]=3<<2;
    size_t i;
    int nBases = 4;
    int nBits = nBases * 2; // => 2^nBits roznych casti, podla vrchnych bitov
    int nParts = myPow(2, nBits);
    int nmask = 0;
    for (i = 0; i < nBits; i++) { // nastavenie "nulovacej" masky -> 00000....1111...1
        nmask |= (1 << i);
    }

    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int KMER_LENGTH = 31;
//    char *input_filename = "/Users/larry/Desktop/tests/raw_part_aa_withoutNs_10";
    char *input_filename = "/Users/larry/Desktop/tests/raw_part_aa";
    char *output_path = "/Users/larry/Desktop/out/"; // TODO: nastavit na input path ak nebude zadana
    char *output_file_name = "output";

    MPI_File m_f;
    MPI_File_open(MPI_COMM_WORLD, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &m_f);
    
    size_t process_part_read_start, process_part_read_end;
    getPartsStartEnd(myrank, npes, input_filename, &process_part_read_start, &process_part_read_end);
//    printf("[%d]: process_part_read_start = %lu => process_part_read_end = %lu\n", myrank, process_part_read_start, process_part_read_end);

    // o rozdeleni informovat vsetky procesy
    size_t _parts_read_start[npes], _parts_read_end[npes];
    size_t *parts_read_start = (size_t*)malloc(sizeof(size_t) * npes);
    size_t *parts_read_end = (size_t*)malloc(sizeof(size_t) * npes);
    for (i = 0; i < npes; i++) {
        parts_read_start[i] = 0;
        parts_read_end[i] = 0;
        
        if (i == myrank) {
            _parts_read_start[i] = process_part_read_start;
            _parts_read_end[i] = process_part_read_end;
        } else {
            _parts_read_start[i] = 0;
            _parts_read_end[i] = 0;
        }
    }
    // aby kazdy vedel zaciatky a konce jednotlivych procesov (treba pre kmers)
    MPI_Allreduce(&_parts_read_start, parts_read_start, npes, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&_parts_read_end, parts_read_end, npes, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
//
//    for(i = 0; i < npes; i++) {
//        printf("[%d]: %lu | %lu -> %lu\n", myrank, i, parts_read_start[i], parts_read_end[i]);
//    }
//    MPI_Finalize();
//    exit(0);
    
    // NMERS COUNTING
    //nmers - sa pocitaju tak ze kazdy process prehladava CAST suboru
    gettimeofday(&start, NULL);
    int *nmers = (int*)malloc(sizeof(int) * nParts);
    for (i = 0; i < nParts; i++) { nmers[i] = 0; }
    getPartNmers(myrank, m_f, process_part_read_start, process_part_read_end, KMER_LENGTH, nBases, nmers);
    gettimeofday(&end, NULL);
    time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
    printf("[%d]: nmers counting time = %lf\n", myrank, time_spent);
    int *totalNmers = (int*)malloc(sizeof(int) * nParts);
    for (i = 0; i < nParts; i++) { totalNmers[i] = 0; }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(nmers, totalNmers, nParts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
//    if (myrank == 0) {
//        for (i = 0; i < nParts; i++) {
//            printf("%lu => %d\n", i, totalNmers[i]);
//        }
//    }
//    MPI_Finalize();
//    exit(0);
    
    // PARTS DETERMINING
    PartLimit *partLimits = (PartLimit*)malloc(sizeof(PartLimit) * nParts);
//    int parts_count_max = npes; // podla poctu procesorov (uzlov, ...)
    int num_parts = 0; // celkovy pocet malych casti (takych aby sa zmestili na jednotlive procesy naraz do pamate)
    countPartsLimits(myrank, npes, nParts, totalNmers, partLimits, &num_parts);
    
//    if (myrank == 0) {
//        for (i = 0; i < num_parts; i++) {
//            printf("part %lu: %d -> %d (%llu)\n", i, partLimits[i].start, partLimits[i].end, partLimits[i].size);
//        }
//    }
//    
    // kmers sa robia tak ze kazdy process prechadza CELY subor (a vybera si len tie pre neho dobre)
    for (i = 0; i < num_parts; i++) {
        int part_id = (int)i;
        int handling_rank = part_id % npes;
        
        if (myrank == handling_rank) {
//            printf("[%d]: handling part %d\n", myrank, part_id);
            // KMERS GETTING
            int *kmers_count = (int*)malloc(sizeof(int));
            *kmers_count = 0;
            uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * partLimits[part_id].count + 100000000); // + rezerva
            gettimeofday(&start, NULL);
            //            printf("[%d]: part_id = %d, expected count = %lu\n", myrank, part_id, partLimits[part_id].count);
            getPartKmers(myrank, m_f, partLimits[part_id].start, partLimits[part_id].end, KMER_LENGTH, nBases, kmers, kmers_count, parts_read_start, parts_read_end);
            //            printf("[%d]: part_id = %d, kmers_count = %d\n", myrank, part_id, *kmers_count);
            gettimeofday(&end, NULL);
            time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: part_id = %d, kmering time = %lf\n", myrank, part_id, time_spent);
            // tu uz mam k-tice z casti
            
            // SORTING
            gettimeofday(&start, NULL);
            uint64_t *aux = (uint64_t*)malloc(sizeof(uint64_t) * (*kmers_count));
            bucketSortSerial(aux, kmers, partLimits[part_id].count, 0, maxLgNBs, 1); // v aux uz budu sortnute? - ne v kmers!!!! (nBits
            free(aux);
            gettimeofday(&end, NULL);
            time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: part_id = %d, sorting time = %lf\n", myrank, part_id, time_spent);
            
//            if (myrank == 0 && part_id == 0) {
//                for (int j = 0; j < 100; j++) {
//                    printf("%llu\n", kmers[j]);
//                }
//            }
            // COUNTING
            gettimeofday(&start, NULL);
            
            char *file_output_counted = (char*) malloc(sizeof(char) * 255);
            sprintf(file_output_counted, "%s%s_counted_%d.bin", output_path, output_file_name, part_id);
            FILE *fout_counted = fopen(file_output_counted, "wb");

            KeyVal *results = (KeyVal*) malloc(sizeof(KeyVal) * (partLimits[part_id].count - (partLimits[part_id].count / 10))); // odhad ze o nieco menej ich bude (unikatnych, spocitanych)
            uint64_t kmerCurr = kmers[0];
            int countCurr = 1;
            int countDist = 0;
            int q;
            for (q = 1; q < partLimits[part_id].count; q++){
                
                if (kmers[q] != kmerCurr) {
                    results[countDist].kmer = kmerCurr;
                    results[countDist].count = countCurr;
                    countDist++;
                    kmerCurr = kmers[q];
                    countCurr = 1;
                } else {
                    countCurr++;
                }
            }
            
            fwrite(results, sizeof(KeyVal), countDist, fout_counted);
            fclose(fout_counted);
            gettimeofday(&end, NULL);
            time_spent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: part_id = %d, writing time = %lf\n", myrank, part_id, time_spent);

            free(kmers);
            free(kmers_count);
        }
    }
    MPI_File_close(&m_f);
    MPI_Finalize();
    
    gettimeofday(&end_total, NULL);
    time_spent = (end_total.tv_sec - start_total.tv_sec) * 1000 + (end_total.tv_usec - start_total.tv_usec)/1000.0;
    printf("[%d]: total time = %lf\n", myrank, time_spent);

    return 0;
}
