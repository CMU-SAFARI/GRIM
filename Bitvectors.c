/*
Jeremie Kim
Modifications to the FastHASH infrastructure to include custom filter.

This filter creates bitvectors for each permutation representing the existence 
of that given permutation in a bin. We split up the genome into a number 
of bins that gives a statistically significant sparsity of sequence existence
across the bins. 

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Bitvectors.h"
#include "Common.h"
#include "RefGenome.h"
/*************************************************/
FILE *_bv_fp = NULL;
FILE *_fbv_fp = NULL;
FILE *_rst_fp = NULL;
unsigned *_bitvector = NULL;
unsigned *_fbitvector = NULL;
char *read_sequence_count = NULL;
char     *_bv_refGenName = NULL;
int      _bv_refGenOff = 0;
unsigned ints_per_perm = 0;
unsigned generate_BV_threshold = 0;
unsigned nuc_threshold = 0;
unsigned SEQUENCE_LENGTH = 100;
double bv_totalCheckTime = 0;
unsigned bitvectorsSize = 0;
long long fbitvectorSize = 0;
unsigned CUR_NUM_BINS = 0;
unsigned iter = 0;
unsigned checkBitvectors_iter = 1;
unsigned last_bin_num = 0;
double avg_checkTime = 0;
long long total_bitvectorsSize = 0;
/*************************************************/

// we want to generate the bitvector file using a reference genome file
void generateIBitvectors(char *refGenome_filename, char *bv_fileName) {
    double startTime = getTime();
    char *refGenName;
    char *refGen;
    int refGenOff = 0;
    //unsigned int bitvectorsSize = 0;
    int flag, l, hv, tmp;
    unsigned int i;
    
    ints_per_perm = pow(4, BV_TOKEN_SIZE)/32;

    _bv_fp = fileOpen(bv_fileName, "w");

    // write the subsequence size used to generate the bitvectors.
    tmp = fwrite(&BV_TOKEN_SIZE, sizeof(BV_TOKEN_SIZE), 1, _bv_fp);
    tmp = fwrite(&BV_BIN_SIZE, sizeof(BV_BIN_SIZE), 1, _bv_fp);
    tmp = fwrite(&BV_NUM_BINS, sizeof(BV_NUM_BINS), 1, _bv_fp);
    tmp = fwrite(&BV_MULTIPLICITY, sizeof(BV_MULTIPLICITY), 1, _bv_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while initializing the bitvectors.\n");
    }
    fprintf(stderr, "Generating Bitvectors from %s", refGenome_filename);
    fflush(stderr);

    //printf("refGenome_filename: %s\n", refGenome_filename);
    if (!initLoadingRefGenome(refGenome_filename)) {
        return;
    }

    char *prev = getMem (CONTIG_NAME_SIZE);
    prev[0]='\0';
    
    do 
    {
        flag = loadRefGenome(&refGen, &refGenName, &refGenOff);
        if (strcmp(prev, refGenName) != 0) {
            fprintf(stderr, "\n - %s ", refGenName);
            fflush(stderr);
            sprintf(prev, "%s", refGenName);
        }
        else {
            fprintf(stderr, ".");
            fflush(stderr);
        }
        
        l = strlen(refGen) - BV_TOKEN_SIZE;

        CUR_NUM_BINS = strlen(refGen) / BV_BIN_SIZE;
        bitvectorsSize = 4*(pow(4, BV_TOKEN_SIZE)/32) * CUR_NUM_BINS * BV_MULTIPLICITY;
        _bitvector = getMem(bitvectorsSize);
        
        for (i = 0; i < l; i++) {
            hv = bv_hashVal(refGen + i);
            if (hv != -1) {
                update_BV(hv, i);
            }
        }

        short len = strlen(refGenName);
        tmp = fwrite(&len, sizeof(len), 1, _bv_fp);
        tmp = fwrite(refGenName, sizeof(char), len, _bv_fp);
        tmp = fwrite(&refGenOff, sizeof(refGenOff), 1, _bv_fp);
        unsigned int refGenLength = strlen(refGen);
        tmp = fwrite(&refGenLength, sizeof(refGenLength), 1, _bv_fp);
        //tmp = fwrite(refGen, sizeof(char), refGenLength, _bv_fp);

        tmp = fwrite(&bitvectorsSize, sizeof(unsigned int), 1, _bv_fp);

        for (i = 0; i < bitvectorsSize/4; i++) {
             tmp = fwrite(&_bitvector[i], sizeof(unsigned int), 1, _bv_fp);
        }
        
        if (tmp == 0) {
            fprintf(stderr, "write error while saving bitvectors,\n");
        }

        // free the bitvector
        freeMem(_bitvector, bitvectorsSize);
        iter++;
    } while (flag);

    finalizeLoadingRefGenome();
    fclose(_bv_fp);

    printf("\nDONE in %0.5fs!\n", (getTime() - startTime));
}

int bv_hashVal(char *seq) {
    int i = 0;
    int val = 0, numericVal = 0;

    while (i < BV_TOKEN_SIZE) {
        //printf("%c", seq[i]);
        switch (seq[i]) {
            case 'A':
                numericVal = 0;
                break;
            case 'C':
                numericVal = 1;
                break;
            case 'G' :
                numericVal = 2;
                break;
            case 'T':
                numericVal = 3;
                break;
            default:
                return -1;
                break;
        }
        val = (val << 2) | numericVal;
        i++;
    }
    //printf("\n");
    return val;
}

// specific to 100 length sequences. 
int update_BV(int hv, int offset) {
    int i;
    int start_bin = (offset - (SEQUENCE_LENGTH + generate_BV_threshold - 1)) / BV_BIN_SIZE;
    if (offset < (SEQUENCE_LENGTH - 1)) { start_bin = 0; }
    int end_bin = offset / BV_BIN_SIZE;
    unsigned BV_offset = hv / ints_per_perm;
    unsigned int_offset = (ints_per_perm - 1) - (hv % ints_per_perm);

    for (i = start_bin; i <= end_bin; i++) {
        unsigned tmp; 
        if (i > (CUR_NUM_BINS - 1)) {
            tmp = ints_per_perm * (CUR_NUM_BINS - 1) + BV_offset;
        }
        else {
            tmp = ints_per_perm * i + BV_offset;
        }
        int j;
        // for multiplicity, we multiply the size of the bitvectors by the multiplicity 
        // and replicate the format. 
        for (j = 0; j < BV_MULTIPLICITY; j++) {
            //tmp += ints_per_perm * (CUR_NUM_BINS);
            int tmp2 = tmp + j * ints_per_perm * CUR_NUM_BINS;
            if (!(1 & (_bitvector[tmp2] >> int_offset))) {
                _bitvector[tmp2] |= (1 << int_offset);
                break;
            }
        }
    }
    return 0;
}

void configBV()
{
    generateBitvectors = &generateIBitvectors;
    loadBitvectors = &loadIBitvectors;
    loadFinalBitvector = &loadIFinalBitvector;
    /*  
        generateHashTable = &generateIHashTable;
        loadHashTable = &loadIHashTable;
        finalizeLoadingHashTable = &finalizeLoadingIHashTable;
        getCandidates = &getIHashTableCandidates;
    */
}

int initLoadingBitvectors(char *fileName) {
    int tmp;
    ints_per_perm = pow(4, BV_TOKEN_SIZE)/32;
    generate_BV_threshold = (SEQUENCE_LENGTH * 5) / 100; // max error threshold of 5%
    nuc_threshold = (SEQUENCE_LENGTH * errThreshold) / 100;
    read_sequence_count = getMem(pow(4, BV_TOKEN_SIZE));
    printf("generate_BV_threshold: %d\n", generate_BV_threshold);
    printf("nuc_threshold: %d\n", nuc_threshold);
    _bv_fp = fileOpen(fileName, "r");
    if (_bv_fp == NULL) {
        return 0;
    }

    tmp = fread(&BV_TOKEN_SIZE, sizeof(BV_TOKEN_SIZE), 1, _bv_fp);
    tmp = fread(&BV_BIN_SIZE, sizeof(BV_BIN_SIZE), 1, _bv_fp);
    tmp = fread(&BV_NUM_BINS, sizeof(BV_NUM_BINS), 1, _bv_fp);
    tmp = fread(&BV_MULTIPLICITY, sizeof(BV_MULTIPLICITY), 1, _bv_fp);
    if (tmp == 0) {
        fprintf(stderr, "read error while getting BV info from BV File\n");
        fflush(stderr);
    }

    printf("BV_TOKEN_SIZE: %d\n", BV_TOKEN_SIZE);
    printf("BV_BIN_SIZE: %d\n", BV_BIN_SIZE);
    printf("BV_NUM_BINS: %d\n", BV_NUM_BINS);

    configBV();
    
    return 1;
}

int initLoadingFinalBitvector(char *fileName) {
    _fbv_fp = fileOpen(fileName, "r");
    if (_fbv_fp == NULL) {
        return 0;
    }
    configBV();
    return 1;
}

int initGeneratingFinalBitvector(char *fileName) {
    printf("%s\n", fileName);
    _fbv_fp = fileOpen(fileName, "w");
    if (_fbv_fp == NULL) {
        return 0;
    }
    return 1;
}

int initGeneratingReadSeqRuntimes(char *fileName) {
    printf("%s\n", fileName);
    _rst_fp = fileOpen(fileName, "w");
    if (_rst_fp == NULL) {
        return 0;
    }
    return 1;
}

void writeRST_listSize (int size) {
    printf("writeRST_listSize size: %d\n", size);
    int tmp = fwrite(&size, sizeof(size), 1, _rst_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while writing the read sequence list size. \n");
    }
}

void writeReadSeqRuntime(double time) {
    int tmp = fwrite(&time, sizeof(time), 1, _rst_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while writing the read sequence runtime. \n");
    }
}

// this is just an array of the bitvector results in order
int loadIFinalBitvector(double *loadTime) {
    double startTime = getTime();
    int tmp = fread(&fbitvectorSize, sizeof(fbitvectorSize), 1, _fbv_fp);
    tmp = fread(_fbitvector, sizeof(unsigned int), (fbitvectorSize+31)/32, _fbv_fp);
    
    printf(" read fbitvectorSize: %llu\n", fbitvectorSize);
    if (tmp == 0) {
        fprintf(stderr, "had problems loading the bitvectors\n");
        return 0;
    }
    
    *loadTime = getTime() - startTime;
    return 1;
}

int loadIBitvectors(double *loadTime) {
    double startTime = getTime();
    unsigned int refGenLength;
    //unsigned int bitvectorsSize;
    int tmp;
    
    short len;
    tmp = fread(&len, sizeof(len), 1, _bv_fp);
    _bv_refGenName = getMem(sizeof(char) * (len + 1));
    tmp = fread(_bv_refGenName, sizeof(char), len, _bv_fp);
    _bv_refGenName[len] = '\0';
    //printf("refGenName: %s\n", _bv_refGenName);

    tmp = fread(&_bv_refGenOff, sizeof(_bv_refGenOff), 1, _bv_fp);
    tmp = fread(&refGenLength, sizeof(refGenLength), 1, _bv_fp);
    //printf("refGenOff: %d\n", _bv_refGenOff);
    //printf("refGenLength: %d\n", refGenLength);
    CUR_NUM_BINS = refGenLength / BV_BIN_SIZE;

    tmp = fread(&bitvectorsSize, sizeof(bitvectorsSize), 1, _bv_fp);
    total_bitvectorsSize += bitvectorsSize;
    _bitvector = getMem(bitvectorsSize);
    tmp = fread(_bitvector, sizeof(unsigned int), bitvectorsSize/4, _bv_fp);
    //printf("bitvectorsSize: %d\n", bitvectorsSize);

    if (tmp == 0) {
        fprintf(stderr, "had problems loading the bitvectors\n");
        return 0;
    }
    
    *loadTime = getTime() - startTime;
    return 1;
}

unsigned* getBitvectors() {
    return _bitvector;
}

// returns 0 if it does pass the threshold and we have to check it DP.
// returns 1 if it does not pass the threshold. 
// returns 2 if if we skip the verification because we already checked the bin. 
unsigned checkBitvectors(unsigned genLoc, char **seq, unsigned subseq_num) {
    memset(read_sequence_count, 0, pow(4, BV_TOKEN_SIZE));
    
    // calculate bin number from genLoc
    unsigned bin_num = (genLoc - (WINDOW_SIZE * subseq_num + nuc_threshold)) / BV_BIN_SIZE;
    // in case of the last bin. (thats longer than the others) 
    if (bin_num >= CUR_NUM_BINS) {
        bin_num = CUR_NUM_BINS - 1;
    }

    if (bin_num == last_bin_num) {
        return 2; // this means we've already seen this happen
    }

    double startTime = 0;
    if (!(checkBitvectors_iter % 1000)) {
        startTime = getTime();
    }

    last_bin_num = bin_num; 
    
    unsigned threshold = SEQ_LENGTH - (nuc_threshold*BV_TOKEN_SIZE);
    // run the filter on that bin against the entire read sequence. 
    unsigned count = 0;
    if (errThreshold == 0) {
        count = 1;
    }
    int i;
    for (i = 0; i < SEQ_LENGTH - BV_TOKEN_SIZE + 1; i++) {
        int hv = bv_hashVal(*seq + i);
        if (hv == -1) return 0;

        int tmp = 0;
        if (BV_MULTIPLICITY != 1) {
            tmp = read_sequence_count[hv]++;
            if (tmp >= BV_MULTIPLICITY) {
                tmp = BV_MULTIPLICITY - 1;
            }
        }

        int tmp_BV_Val = _bitvector[(tmp * (CUR_NUM_BINS) * ints_per_perm) + 
                                    (bin_num * ints_per_perm) + (hv / ints_per_perm)]
                                    >> ((ints_per_perm - 1) - (hv % ints_per_perm));

        if (errThreshold == 0) {
            count &= (1 & tmp_BV_Val);
            if (count == 0) {
                break;
            }
        }
        else {
            count += (1 & tmp_BV_Val);
            if (count >= threshold) {
                break;
            }
        }
    }
    
    if (!(checkBitvectors_iter % 1000)) {
        avg_checkTime += (getTime() - startTime);
    }
    checkBitvectors_iter++;

    if (errThreshold == 0) {
        return !count;
    }
    else {
        return count < threshold;
    }
}

double get_bv_totalCheckTime() {
    return avg_checkTime * 1000; //bv_totalCheckTime;
}

long long get_total_bitvectorsSize() {
    return total_bitvectorsSize;
}

/*
void clear_fbv() {
    
}

// this will write the final bitvector to the file. 
int write_fbv() {

}


int initialize_fbv(unsigned NUM_BINS, unsigned num_reads) {
    unsigned size = ((NUM_BINS + 7) / 8) * num_reads;
}



char get_fbv_result(unsigned fbv_index) {
    if (fbv_index >= fbitvectorSize) {
        return 0;   
    }

    return 1 & (_fbitvector[fbv_index / 32] >> (31 - (fbv_index % 31))); 
}

// writing space that will eventually be overwritten
void fbv_write_size_space() {
    long long tmp2 = 0;
    int tmp = fwrite(&tmp2, sizeof(tmp2), 1, _fbv_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while writing data to fbv\n");
    }
}

void fbv_write_size(long long size) {
    fseek(_fbv_fp, -(((size + 7)/8) + (sizeof(long long))), SEEK_CUR);
    int tmp = fwrite(&size, sizeof(size), 1, _fbv_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while writing the size of the fbv\n");
    }
    printf(" write fbitvectorSize: %llu\n", size);
    fseek(_fbv_fp, 0, SEEK_END);
}

void fbv_write_data(unsigned data) {
    int tmp = fwrite(&data, sizeof(data), 1, _fbv_fp);
    if (tmp == 0) {
        fprintf(stderr, "Write error while writing data to fbv\n");
    }
}
*/
