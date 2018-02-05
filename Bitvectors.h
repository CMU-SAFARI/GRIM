#ifndef __BITVECTORS__
#define __BITVECTORS__

void generateIBitvectors(char *refGenome_filename, char *bv_filename);
int  bv_hashVal(char *seq);
int  update_BV(int hv, int offset);
void configBV();
int initLoadingBitvectors(char *fileName);
int initLoadingFinalBitvector(char *fileName);
int initGeneratingFinalBitvector(char *fileName);
int initGeneratingReadSeqRuntimes(char *fileName);
void writeRST_listSize(int size);
void writeReadSeqRuntime(double time);
int loadIBitvectors(double *loadTime);
int loadIFinalBitvector(double *loadTime);
unsigned* getBitvectors();
unsigned checkBitvectors(unsigned genLoc, char** seq, unsigned subseq_num);
double get_bv_totalCheckTime();
long long get_total_bitvectorsSize();
char get_fbv_result(unsigned fbv_index);
void fbv_write_size_space();
void fbv_write_size(long long size);
void fbv_write_data(unsigned data);

void (*generateBitvectors)(char *refGenome_filename, char *bv_filename);
int  (*loadBitvectors)(double *loadTime);
int  (*loadFinalBitvector)(double *loadTime);
#endif
