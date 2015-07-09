#include <mpi.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

// 128MB = 134217728
//
#define MAX_MEMORY_PROCESSOR_FILE_LOAD 1073741824//1073741824// max pamate pre proces/processor - 1 GB (1 073 741 824) - najviac tolkoto nacitat naraz
        // maximalne nieco pod 2GB, lebo INT
#define PART_SIZE_MAX 1073741824//134217728 //1073741824
#define MAX_READS_AMOUNT 100000000 // TODO (ked je malo -> segmentation fault)
#define MAX_READ_SIZE 1000
#define KMERS_BUFFER_MAX 100 // max pocet k-tic kolko sa musi nazbierat aby sa zapisali do suboru

char BASES[255];
char BASES_2[255];

//#define baseBitShift 8 // nBits   // pocet bucketov load balancing 2^n  (BASEBITSHIFT je to uvodne delenie - je tam 8 -> 2exp8 = 256, odporucam 16*pocet jadier)
#define maxLgNBs 12				//max. pocet bucketov v deleni - 2^n (MAXLGNBSS sa tyka bucket sortu, da sa pskusat ci 11, 12, 13, 14 vyssie by som nesiel, ale nastavene je to na 12, ze to takdobre chodilo)

typedef struct PartLimit {
    int start;
    int end;
    uint64_t num_parts;
    uint64_t size;
    uint64_t count;
} PartLimit;

typedef struct KeyVal {
    uint64_t kmer;
    int count;
} KeyVal;


int init(int argc, char *argv[]/*, int *READ_LENGTH, int *LINE_LENGTH,*/, int *KMER_LENGTH/*, int *KMERS_IN_READ*/, int *INPUT_FILE_TYPE, char **input_file_full_path, char **output_path, char **output_file_name) {
    BASES[65]=0; BASES[67]=1; BASES[71]=2; BASES[84]=3;
    BASES_2[65]=0<<2; BASES_2[67]=1<<2; BASES_2[71]=2<<2; BASES_2[84]=3<<2;
    
    // PARAMETRE
    char c;
    int index;
    opterr = 0;
    while ((c = getopt(argc, argv, "k:o:ft")) != -1) {
        switch (c) {
            case 'l':                               // -l READ_LENGTH
//                *READ_LENGTH = atoi(optarg);
                break;
            case 'k':                               // -k KMER_LENGTH
                *KMER_LENGTH = atoi(optarg);
                break;
            case 'o':                               // -o output_path
                *output_path = optarg;
                break;
            case 'f':                               // -f output_filename
                *output_file_name = optarg;
                break;
            case 't':                               // -t filetype (0 = raw, 1 = fasta)
                if (strcmp(optarg, "raw") == 0) {
                    *INPUT_FILE_TYPE = 0;
                } else {
                    *INPUT_FILE_TYPE = 1;
                }
                break;
//            case 'i':
//                printf("%s\n", optarg);
//                *input_file_full_path = optarg;
//                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
            default:
                abort();
                
        }
    }
    
    for (index = optind; index < argc; index++) { // ostatne non-option argumenty
//        printf ("Non-option argument %d = %s\n", index, argv[index]);
        *input_file_full_path = argv[index];
    }
//    *KMERS_IN_READ = *READ_LENGTH - *KMER_LENGTH + 1;
//    *LINE_LENGTH = *READ_LENGTH + 1;
    return 0;
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

// ready nechat v povodnom raw_buffer, zmazat nove riadky, a ziskat indexy (pozicie readov v raw_buffer) - tj. ich dlzky
void getReadsFromRaw(char *raw_buffer, size_t raw_buffer_size, uint64_t *indexes, uint64_t *num_of_reads) {
    int i, indexCurrent = 0, indexCurrentId = 0;

    char *currSeq = (char*)malloc(sizeof(char) * MAX_READ_SIZE);
    int currSeqLength = 0;
    
    for (i = 0; i < raw_buffer_size; i++) {

        if (raw_buffer[i] == '\n' || raw_buffer[i] == EOF || raw_buffer == '\0') { // dalsi read
            strncpy(&raw_buffer[indexCurrent], currSeq, currSeqLength); // prepisat do povodneho buffera (zo zaciatku (indexCurrent))
            indexes[indexCurrentId] = indexCurrent; // zaciatocny index
            indexCurrent += currSeqLength + 1;
            raw_buffer[indexCurrent - 1] = '\0';
            indexCurrentId++;
            (*num_of_reads)++;
            memset(currSeq, '\0', MAX_READ_SIZE);
            currSeqLength = 0;
            continue;
        } else {
            currSeq[currSeqLength] = raw_buffer[i];
            currSeqLength++;
        }
    }
//    strncpy(&raw_buffer[indexCurrent], currSeq, currSeqLength); // prepisat do povodneho buffera (zo zaciatku (indexCurrent))
    indexes[indexCurrentId] = indexCurrent; // (aby sme mali aj posledny index)
//    indexCurrent += currSeqLength + 1;
//    raw_buffer[indexCurrent - 1] = '\0';
//    indexCurrentId++;
//    (*num_of_reads)++;
//    memset(currSeq, '\0', MAX_READ_SIZE);
//    currSeqLength = 0;
}


// ready necham v povodnom 'fasta_buffer', kde zmazem akurat headery a ostatne info (necham len sekvencie) a do indexes dam pozicie zaciatkov readov
void getReadsFromFasta(char *fasta_buffer, uint64_t fasta_buffer_size, uint64_t *indexes, uint64_t *num_of_reads) {
    uint64_t i, indexCurrent = 0;
    uint64_t indexCurrentId = 0;
    int isHeaderLine = 1;
    int isSequenceLine = 0;
    
    char *currSeq = (char*)malloc(sizeof(char) * MAX_READ_SIZE);
    int currSeqLength = 0;
    
    for (i = 0; i < fasta_buffer_size; i++) {
        if (fasta_buffer[i] == '\0')
            break;
        //        printf("%d(%c)\n", i, fasta_buffer[i]);
        if (fasta_buffer[i] == '>' || fasta_buffer[i] == '@' || fasta_buffer[i] == EOF) {
            
            if (isSequenceLine) {
                //                printf("%d: %s\n\n", currSeqLength, currSeq);
                strncpy(&fasta_buffer[indexCurrent], currSeq, currSeqLength); // prepisat do povodneho buffera (zo zaciatku (indexCurrent))
                indexes[indexCurrentId] = indexCurrent; // zaciatocny index
                indexCurrent += currSeqLength + 1;
                fasta_buffer[indexCurrent - 1] = '\0';
                indexCurrentId++;
                (*num_of_reads)++;
                memset(currSeq, '\0', MAX_READ_SIZE);
                currSeqLength = 0;
            }
            isHeaderLine = 1;
            isSequenceLine = 0;
            continue;
        }
        
        if (fasta_buffer[i] == '+') { // vynechat kvlalitu
            isHeaderLine = 0;
            isSequenceLine = 0;
            continue;
        }
        
        if (isHeaderLine && fasta_buffer[i] == '\n') { // po headeri ide sekvencia
            isHeaderLine = 0;
            isSequenceLine = 1;
            continue;
        }
        
        if (isSequenceLine) {
            
            if (fasta_buffer[i] != '\n') {
                currSeq[currSeqLength] = fasta_buffer[i];
                currSeqLength++;
                continue;
            }
        }
    }
    //    printf("%d: %s\n\n", currSeqLength, currSeq);
    strncpy(&fasta_buffer[indexCurrent], currSeq, currSeqLength);
    indexes[indexCurrentId] = indexCurrent; // zaciatocny index
//    printf("lastIndex = %llu, value = %llu\n", indexCurrentId, indexes[indexCurrentId]);
    indexCurrent += currSeqLength + 1;
    fasta_buffer[indexCurrent - 1] = '\0';
//    indexCurrentId++;
//    (*num_of_reads)++; // TODO: chyba - bude chybat posledny read?
    memset(currSeq, ' ', MAX_READ_SIZE);
    currSeqLength = 0;
}


// najde offset zaciatku posledneho zaznamu(readu) podla typu vstupneho suboru
size_t getRestOffset(int INPUT_FILE_TYPE, char *buf, size_t buf_size) {
    size_t i, offset = 0;
    
    if (INPUT_FILE_TYPE == 0) { // raw
        for (i = buf_size-1; buf[i] != '\n'; i--) { // novy riadok == novy zaznam
            offset++;
        }
        if (i < buf_size-1) { // iba ak sa to vobec skusalo
            //offset++; // pridat jednu poziciu ?
        }
    } else if (INPUT_FILE_TYPE == 1) { // fasta
        for (i = buf_size-1; (buf[i] != '>' && buf[i] != '@'); i--) { // '>', '@' == novy zaznam
            offset++;
        }
    }
    return offset;
}

void getNmers(char *filename, int INPUT_FILE_TYPE, int nBases, int KMER_LENGTH, int *nmers) {
    int nBits = nBases * 2; // => 2^nBits roznych casti, podla vrchnych bitov
    double nPartsDouble = pow(2.0, (double)nBits);
    int nParts = (int)nPartsDouble;
    int nmask = 0;
    uint64_t i;
    
    for (i = 0; i < nParts; i++) { // vynulovanie pola
        nmers[i] = 0;
    }
    for (i = 0; i < nBits; i++) { // nastavenie "nulovacej" masky -> 00000....1111...1
        nmask |= (1 << i);
    }
    FILE *f = fopen(filename, "r");
    
    size_t rest_max = 2048; // maximalny zvysok (ak sa neda subor nacitat naraz a na konci jedneho citania nemame kompletny read)
    size_t rest_offset = 0; // kolko pridat pamate z minuleho kola (zostatok z minuleho citania (posledny read))
    uint64_t count_read = 0;

    char *buf = (char *) malloc(sizeof(char) * MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
    if (buf == NULL) printf("NEPODARILO sa alokovat *buf afioaj\n");
    memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
    
    char *rest = (char*)malloc(sizeof(char) * rest_max); // alokovat pamat pre zvysok (posledny read)
    if (rest == NULL) printf("NEPODARILO sa alokovat *rest fasfgge\n");
    memset(rest, '\0', MAX_READ_SIZE);
    
    uint64_t *indexes = (uint64_t*)malloc(sizeof(uint64_t) * MAX_READS_AMOUNT); //max pocet readov
    if (indexes == NULL) printf("NEPODARILo sa alokovat *indexes gtertiort\n");

    uint64_t *num_of_reads = (uint64_t*)malloc(sizeof(uint64_t));
    *num_of_reads = 0;
    
    char *read = (char*) malloc(sizeof(char) * MAX_READ_SIZE);
    memset(read, '\0', MAX_READ_SIZE);

    count_read = fread(buf, sizeof(char), MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max, f);
    uint64_t nmer_value = -1;

    while (count_read > 0) {
        //printf("count_read = %d\n", count_read);
        if (rest_offset > 0) {
            strncpy(&buf[0], &rest[0], rest_offset);// pridat zvysok na zaciatok buf
            memset(rest, '\0', rest_offset); // vynulovat aktualny zvysok
        }
        size_t old_rest_offset = rest_offset; // potrebujeme stary offset aby sme vedeli kolko sa posunut
        
        rest_offset = getRestOffset(INPUT_FILE_TYPE, buf, count_read + old_rest_offset);
        if (rest_offset > 0) {
            strncpy(rest, &buf[count_read + old_rest_offset - rest_offset], rest_offset); // zapamatat si zvysok
        }
        memset(&buf[count_read + old_rest_offset - rest_offset], '\0', rest_offset);// odstranit zvysok z konca buf
        
        if (INPUT_FILE_TYPE == 0) { // raw
            getReadsFromRaw(buf, count_read + old_rest_offset, indexes, num_of_reads);
        } else {
            getReadsFromFasta(buf, count_read + old_rest_offset, indexes, num_of_reads);
        }
//        printf("count_read + old_rest_offset = %llu\n", count_read + old_rest_offset);
//        printf("%s", buf);
//        printf("i = %llu, num_of_reads = %llu\n", i, *num_of_reads);
//        printf("count_read = %llu\n", count_read);
        int k;
//        printf("buffer_size = %llu\n", sizeof(char) * MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
        for (k = 0; k < (*num_of_reads); k++) {
            int read_length = (int)(indexes[k+1]-1 - indexes[k]);
            memcpy(read, &buf[ indexes[k] ], read_length);
            if (strlen(read) < KMER_LENGTH || strchr(read, 'N')) {
                continue;
            }
            read[read_length] = '\0';
            int first_nmer = 0;
            
            for (i = 0; i < nBases; i++) {
                first_nmer |= (BASES[read[i]] << (nBits - i * 2 - 2));
            }
            nmers[first_nmer]++;
            int nmer = first_nmer;
            for (i = nBases; i < read_length - KMER_LENGTH + nBases; i++) {
                nmer = ((nmer << 2) & nmask) + BASES[read[i]];
                nmers[nmer]++;
            }
            memset(read, '\0', MAX_READ_SIZE);
        }
        nmer_value = -1;
        *num_of_reads = 0;
        memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
        count_read = fread(buf, sizeof(char), MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max, f);
        
    } // END while (count_read > 0)
    
    free(read); // TODO
    
    free(num_of_reads);
    free(indexes);
    free(buf);
    free(rest);
    
    fclose(f);
}

int main(int argc, char *argv[]) {
    const int debug = 0;
    
    struct timeval start, end, start_total, end_total;
    gettimeofday(&start_total, NULL);
    double timeSpent;
    
//    int READ_LENGTH = 76, KMERS_IN_READ = 46, LINE_LENGTH = 77;
    int KMER_LENGTH = 31;
    int INPUT_FILE_TYPE = 1; // 0: raw reads, 1: fasta
    char *input_file_full_path = "";
    char *output_path = ""; // nastavit na input path ak nebude zadana
    char *output_file_name = "output";

    init(argc, argv/*, &READ_LENGTH, &LINE_LENGTH*/, &KMER_LENGTH/*, &KMERS_IN_READ*/, &INPUT_FILE_TYPE, &input_file_full_path, &output_path, &output_file_name);
    
    int npes, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int PARTS = npes;

    if (myrank == 0) {
//        printf("READ_LENGTH = %d\n", READ_LENGTH);
//        printf("LINE_LENGTH = %d\n", LINE_LENGTH);
        printf("KMER_LENGTH = %d\n", KMER_LENGTH);
//        printf("KMERS_IN_READ = %d\n", KMERS_IN_READ);
        printf("INPUT_FILE_TYPE = %d\n", INPUT_FILE_TYPE);
        
        printf("input_file_full_path = %s\n", input_file_full_path);
        printf("output_path = %s\n", output_path);
        printf("output_file_name = %s\n", output_file_name);
    }
    FILE *f_test = fopen(input_file_full_path, "r");
    if (f_test == NULL) {

        if (myrank == 0) {
            fprintf(stderr, "Error opening file: '%s'\n", input_file_full_path);
        }
        MPI_Finalize();
        exit(0);
    }
    fclose(f_test);
    
    uint64_t i;

    MPI_Status m_status;
    // nmers -> vrchne bity (k vypoctu ich frekvencii)
    int nBases = 4;
    int nBits = nBases * 2; // => 2^nBits roznych casti, podla vrchnych bitov

    double nPartsDouble = pow(2.0, (double)nBits);
    int nParts = (int)nPartsDouble;//;pow(2.0, (double)nBits);
    
    int *nmers = (int*) malloc(sizeof(int) * nParts);
    int nmask = 0;
    for (i = 0; i < nParts; i++) { // vynulovanie pola
        nmers[i] = 0;
    }
    for (i = 0; i < nBits; i++) { // nastavenie "nulovacej" masky -> 00000....1111...1
        nmask |= (1 << i);
    }
    
    uint64_t KMER_SIZE = sizeof(uint64_t); // 8B
    
    if (myrank == 0) {
        gettimeofday(&start, NULL);
        getNmers(input_file_full_path, INPUT_FILE_TYPE, nBases, KMER_LENGTH, nmers);
        gettimeofday(&end, NULL);
        timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
        printf("nmers counting: %lf\n", timeSpent);
        
        for (i = 1; i < npes; i++) { // poslat kazdemu spocitane nmery
            MPI_Send(nmers, nParts, MPI_INT, (int)i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(nmers, nParts, MPI_INT, 0, 0, MPI_COMM_WORLD, &m_status); // ostatni prijmu tuto informaciu
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
//    MPI_Finalize();
//    exit(0);
//    
    MPI_Offset m_filesize;
    MPI_File m_file;
    MPI_File_open(MPI_COMM_WORLD, input_file_full_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &m_file);
    MPI_File_get_size(m_file, &m_filesize);
    
    
    // ZISTOVANIE JEDNOTLIVYCH CASTI (parts)
    // urcite ktore k-tice pojdu kam
    //    PartLimit partLimits[nParts];
    gettimeofday(&start, NULL);
    PartLimit *partLimits = (PartLimit*)malloc(sizeof(PartLimit) * nParts);
    //    printf("sizeof(PartLimit) = %d\nsizeof(partLimits) = %d\n", sizeof(PartLimit), sizeof(PartLimit) * nParts);
    int parts_count_max = npes; // podla poctu procesorov (uzlov, ...)
    
    // bolo by dobre aby casti boli rovnako velke, casti je 'parts_count_max'
    uint64_t kmers_count_total = 0;
    for (i = 0; i < nParts; i++) {
        kmers_count_total += nmers[i];
    }
    //    printf("kmers_count_total = %llu\n", kmers_count_total);
    // rozdelit na rovnomerne casti, ked sa nezmesti do pamate (vacsie ako part_size_max) tak na viac kol
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
    int num_parts = 0;
    int num_parts_current = 0;
    int part_size_curr = 0;
    int part_count_curr = 0;
    partLimits[num_parts].start = 0;
    partLimits[num_parts].end = 0;
    for (i = 0; i < nParts; i++) {
        part_count_curr += nmers[i];
        part_size_curr += nmers[i] * KMER_SIZE;
        num_parts_current++;
        
        if (part_size_curr > part_size_equal || i == nParts - 1) {
            partLimits[num_parts].end = (int)i;
            partLimits[num_parts].size = part_size_curr;
            partLimits[num_parts].count = part_count_curr;
            partLimits[num_parts].num_parts = num_parts_current;
            //            printf("[%d] set id = %d, start = %d, end = %d\n", myrank, num_parts, partLimits[num_parts].start, partLimits[num_parts].end);
            num_parts++;
            partLimits[num_parts].start = (int)i + 1;
            partLimits[num_parts].end = (int)i + 1;
            part_size_curr = 0;
            num_parts_current = 0;
            part_count_curr = 0;
        }
    }
    gettimeofday(&end, NULL);
    timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
//    printf("[%d]: parts counting: %lf\n", myrank, timeSpent);
    gettimeofday(&start, NULL);

 
    if (debug && myrank == 0) {
        printf("number of parts = %d\n", num_parts);
        for (i = 0; i < num_parts; i++) {
            printf("part %d: %d -> %d (%llu)\n", (int)i, partLimits[i].start, partLimits[i].end, partLimits[i].size);
        }
    }
//    MPI_Finalize();
//    exit(0);
    // TTTTTTTTTTTTTTT
    // UUUUUUUUUUUUUUU
    // ZAPIS JEDNOTLIVYCH CASTI (parts) DO ROZNYCH SUBOROV
//    printf("[%d]: END ROUND %d, start_level = %d, end_level = %d, count = %llu, size = %llu, timeSpent = %lf\n", myrank, round, start_level, end_level, partLimits[part_id].count, partLimits[part_id].size, timeSpent);

    for (i = 0; i < num_parts; i++) {
        int handling_rank = (int)i % npes;

        if (myrank == handling_rank) {
//            printf("[%d]: handling rank = %d, handling part = %d\n", myrank, handling_rank, i);
            int start_level = (partLimits[i]).start;
            int end_level = (partLimits[i]).end;
            printf("[%d]: part %d: start_level = %d -> end_level = %d (size = %llu)\n", myrank, (int)i, start_level, end_level, partLimits[i].size);
            
            char *file_output = (char*) malloc(sizeof(char) * 255);
            sprintf(file_output, "%stmp_%s_%d.bin", output_path, output_file_name, (int)i);
            FILE *fout = fopen(file_output, "wb");
           
            int rest_max = 2048; // maximalny zvysok (ak sa neda subor nacitat naraz a na konci jedneho citania nemame kompletny read)
            char *buf = (char *) malloc(sizeof(char) * MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
            if (buf == NULL) printf("NEPODARILO sa alokovat *buf afioaj\n");
            memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
            
            char *rest = (char*)malloc(sizeof(char) * rest_max); // alokovat pamat pre zvysok (posledny read)
            if (rest == NULL) printf("NEPODARILO sa alokovat *rest fasfgge\n");
            memset(rest, '\0', MAX_READ_SIZE);
            uint64_t rest_offset = 0; // kolko pridat pamate z minuleho kola (zostatok z minuleho citania (posledny read))
            int count_read = 0;
            
            MPI_File_seek(m_file, 0, MPI_SEEK_SET); // nastavit na zaciatok
            MPI_File_read(m_file, &buf[rest_offset], MAX_MEMORY_PROCESSOR_FILE_LOAD, MPI_CHAR, &m_status);
            MPI_Get_count(&m_status, MPI_CHAR, &count_read);
            
            int k, q = 0;
            uint64_t *indexes = (uint64_t*)malloc(sizeof(uint64_t) * MAX_READS_AMOUNT); //max pocet readov
            if (indexes == NULL) printf("NEPODARILo sa alokovat *indexes gtert iort\n");
            uint64_t *num_of_reads = (uint64_t*)malloc(sizeof(uint64_t));
            *num_of_reads = 0;
            
            int kmers_buffer_max = KMERS_BUFFER_MAX; // TODO
            uint64_t *kmers_buffer = (uint64_t*) malloc(sizeof(uint64_t) * (kmers_buffer_max + 1));
            uint64_t kmers_counter_current = 0;
            
            char *read = (char*) malloc(sizeof(char) * MAX_READ_SIZE);
            memset(read, '\0', MAX_READ_SIZE);

            int nmer_value = -1;
            uint64_t kmer = 0;
            int ii;
            while (count_read > 0) {
                printf("[%d]: count_read = %d\n", myrank, count_read);
                if (rest_offset > 0) {
                    strncpy(&buf[0], &rest[0], rest_offset);// pridat zvysok na zaciatok buf
                    memset(rest, '\0', rest_offset); // vynulovat aktualny zvysok
                }
                uint64_t old_rest_offset = rest_offset; // potrebujeme stary offset aby sme vedeli kolko sa posunut
                rest_offset = getRestOffset(INPUT_FILE_TYPE, buf, count_read + old_rest_offset);
                if (rest_offset > 0) {
                    strncpy(rest, &buf[count_read + old_rest_offset - rest_offset], rest_offset); // zapamatat si zvysok
                }
                memset(&buf[count_read + old_rest_offset - rest_offset], '\0', rest_offset);// odstranit zvysok z konca buf
                
                if (INPUT_FILE_TYPE == 0) { // raw
                    getReadsFromRaw(buf, count_read + old_rest_offset, indexes, num_of_reads);
                } else if (INPUT_FILE_TYPE == 1) {
                    getReadsFromFasta(buf, count_read + old_rest_offset, indexes, num_of_reads);
                }
                
                //                printf("[%d]: i = %lu, num_of_reads = %d\n", myrank, i, *num_of_reads);
                for (k = 0; k < (*num_of_reads); k++) {
//                    printf("index %d = %d - %d \t:\t%s\n", k, indexes[k], indexes[k+1]-1, &buf[ indexes[k] ]);
//                    continue;
                    int read_length = (int)(indexes[k+1]-1 - indexes[k]);
                    memcpy(read, &buf[ indexes[k] ], read_length);
                    if (strlen(read) < KMER_LENGTH || strchr(read, 'N')) {
                        continue;
                    }
//                    printf("%s\n", read);
                    for (ii = 0; ii < KMER_LENGTH; ii++) {
                        kmer |= (((uint64_t)BASES[read[ii]]) << (62 - ii * 2)); //TODO: check ?
                    }
                    nmer_value = (kmer >> (64 - nBits)); // ziskat hodnotu prvych nBits tejto k-tice

                    if (nmer_value >= start_level && nmer_value <= end_level) {
                        kmers_buffer[kmers_counter_current] = kmer;
                        kmers_counter_current++;
                        
                        if (kmers_counter_current == kmers_buffer_max) {
                            
                            if (!debug) {
                                fwrite(kmers_buffer, KMER_SIZE, kmers_counter_current, fout);
                            }
                            kmers_counter_current = 0;
                        }
                    }
                    for (ii = 0; ii < (read_length - KMER_LENGTH + 1); ii++) { // TODO: na ostatnych miestach treba dat -1 !!!
                        kmer = (uint64_t)((kmer << 2) | ((uint64_t)BASES_2[read[KMER_LENGTH + ii]])); // << 2
                        nmer_value = (kmer >> (64 - nBits)); // ziskat hodnotu prvych nBits tejto k-tice

                        if (nmer_value >= start_level && nmer_value <= end_level) {
                            kmers_buffer[kmers_counter_current] = kmer;
                            kmers_counter_current++;

                            if (kmers_counter_current == kmers_buffer_max) {

                                if (!debug) {
                                    fwrite(kmers_buffer, KMER_SIZE, kmers_counter_current, fout);
                                }
                                kmers_counter_current = 0;
                            }
                        }
                        
                    }
                    memset(read, '\0', MAX_READ_SIZE);
                    nmer_value = -1;
                    kmer = 0;
                    
                } // for k < (*num_of_reads)
                
                if (kmers_counter_current < kmers_buffer_max && kmers_counter_current > 0) { // zapisat zvysok

                    if (!debug) {
                        fwrite(kmers_buffer, KMER_SIZE, kmers_counter_current, fout);
                    }
                    kmers_counter_current = 0;
                }
                *num_of_reads = 0;
                memset(buf, '\0', MAX_MEMORY_PROCESSOR_FILE_LOAD + rest_max);
                MPI_File_read(m_file, &buf[rest_offset], MAX_MEMORY_PROCESSOR_FILE_LOAD, MPI_CHAR, &m_status);
                MPI_Get_count(&m_status, MPI_CHAR, &count_read);
                q++;
                
            } // END while (count_read > 0)

            free(read); // TODO
            free(kmers_buffer); // TODO
            
            free(num_of_reads);
            free(indexes);
            free(buf);
            free(rest);
            
//            MPI_File_close(&m_fout);
            free(file_output);
            fclose(fout);
//            fclose(fout_txt);
        }
    } // END // zapis jednotlivych casti do suborov
    
    gettimeofday(&end, NULL);
    timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;

    if (!debug) {
        printf("[%d]: zapis casti do suborov %lf\n", myrank, timeSpent);
    }

    MPI_File_close(&m_file);
//    MPI_Finalize();
//    exit(0);
    
    /// QQQQQQQ
    // kazda cast je samostatne nacitatelne a zmesti sa do pamati

    for (i = 0; i < num_parts; i++) {
        int part_id = (int)i;
        int handlingRank = part_id % npes;

        
        if (myrank == handlingRank) {
            char *file_input = (char*) malloc(sizeof(char) * 255);
            sprintf(file_input, "%stmp_%s_%d.bin", output_path, output_file_name, part_id); // vytvoreny tmp subor
//            printf("working on %s\n", file_input);
            
            MPI_File_open(MPI_COMM_WORLD, file_input, MPI_MODE_RDONLY, MPI_INFO_NULL, &m_file);
            MPI_File_get_size(m_file, &m_filesize);
            
    //        printf("[%d]: part_id = %d, filesize=%llu\n", myrank, part_id, m_filesize);
    //        printf("[%d]: part_id = %d, part.size = %llu, part.count = %llu\n", myrank, part_id, partLimits[part_id].size, partLimits[part_id].count);
            
            gettimeofday(&start, NULL);
            uint64_t *kmers_all = (uint64_t*) malloc(sizeof(uint64_t) * partLimits[part_id].count);
            uint64_t *aux = (uint64_t*) malloc(sizeof(uint64_t) * partLimits[part_id].count);
            MPI_File_read(m_file, kmers_all, (int)partLimits[part_id].count, MPI_UINT64_T, &m_status);
            MPI_File_close(&m_file); // uz to mame v kmers_all -> mozme zavriet
            gettimeofday(&end, NULL);
            timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: load %llu kmerov: %lf\n", myrank, partLimits[part_id].count, timeSpent);
/*
                // po castiach sortovat? (ALEBO TO SORTNUT CELE NARAZ) ?  - cele naraz sa oplati viac
                int ii;
                gettimeofday(&start, NULL);
                for (ii = partLimits[part_id].start; ii <= partLimits[part_id].end; ii++) {
                    uint64_t *kmers = (uint64_t*) malloc(sizeof(uint64_t) * nmers[ii]);
                    uint64_t *aux = (uint64_t*) malloc(sizeof(uint64_t) * nmers[ii]);
                    
                    int j, k = 0;
                    uint64_t nmer_value;
                    printf("doing: %d\n", ii);
                    for (j = 0; j < partLimits[part_id].count; j++) { // zo vsetkych vybrat prave 'ii' nmery
                        nmer_value = (kmers_all[j] >> (64 - nBits));

                        if ((int)nmer_value == ii) {
                            kmers[k] = kmers_all[j];
                            k++;
                        }
                    }
                    printf("%d done...\n", ii);
                    bucketSortSerial(aux, kmers, nmers[ii], nBits, maxLgNBs, 1); // v aux uz budu sortnute? - ne v kmers!!!!
                    printf("%d sorting done...\n", ii);
                }
*/

            gettimeofday(&start, NULL);
    //        bucketSortSerial(aux, kmers, partLimits[part_id].count, nBits, maxLgNBs, 1);
            bucketSortSerial(aux, kmers_all, partLimits[part_id].count, 0, maxLgNBs, 1); // v aux uz budu sortnute? - ne v kmers!!!! (nBits -> mocniny 2 - posuvat)
            free(aux);
            
            gettimeofday(&end, NULL);
            timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: sorting: %lf\n", myrank, timeSpent);

            // SORTED
//            printf("writing sorted...\n");
//            char *file_output_sorted = (char*) malloc(sizeof(char) * 255);
//            sprintf(file_output_sorted, "%s%s_sorted_%d.bin", output_path, output_file_name, part_id);
//            FILE *fout_sorted = fopen(file_output_sorted, "wb");
//            fwrite(kmers_all, KMER_SIZE, partLimits[part_id].count, fout_sorted);
//            fclose(fout_sorted);
            
            // COUNTING
//            printf("[%d]: part = %d, zapis counted: ...\n", myrank, part_id);
            gettimeofday(&start, NULL);
            FILE *fout_counted;
            if (!debug) {
                char *file_output_counted = (char*) malloc(sizeof(char) * 255);
                sprintf(file_output_counted, "%s%s_counted_%d.bin", output_path, output_file_name, part_id);
                fout_counted = fopen(file_output_counted, "wb");
            }
//            char *file_output_counted_txt = (char*) malloc(sizeof(char) * 255);
//            sprintf(file_output_counted_txt, "%s%s_counted_%d.txt", output_path, output_file_name, part_id);
//            FILE *fout_counted_txt = fopen(file_output_counted_txt, "w");
            
            KeyVal *results = (KeyVal*) malloc(sizeof(KeyVal) * partLimits[part_id].count);
            uint64_t kmerCurr = kmers_all[0];
            int countCurr = 1;
            int countDist = 0;
            int q;
            for (q = 1; q < partLimits[part_id].count; q++){
                
                if (kmers_all[q] != kmerCurr) {
                    results[countDist].kmer = kmerCurr;
                    results[countDist].count = countCurr;
                    countDist++;
                    kmerCurr = kmers_all[q];
                    countCurr = 1;
                } else {
                    countCurr++;
                }
            }

            if (!debug) {
                fwrite(results, sizeof(KeyVal), countDist, fout_counted);
                fclose(fout_counted);
            }
//            fclose(fout_counted_txt);
            gettimeofday(&end, NULL);
            timeSpent = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec)/1000.0;
            printf("[%d]: counting: %lf\n", myrank, timeSpent);
            
            free(kmers_all);
        }
        
    }
    
    gettimeofday(&end_total, NULL);
    double time_elapsed = (end_total.tv_sec - start_total.tv_sec) * 1000 + (end_total.tv_usec - start_total.tv_usec)/1000.0;
    printf("[%d] time elapsed=%lfms\n", myrank, time_elapsed);

    // zmazat tmp subory
    
    
    MPI_Finalize();
//    exit(0);
    return 0;
}
