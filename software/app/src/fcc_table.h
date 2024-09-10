#ifndef FCC_FUNCTIONS_H
#define FCC_FUNCTIONS_H

#define MAX_LINE_LENGTH 256
#define MAX_SERVICE_LENGTH 100
#define MAX_NOTES_LENGTH 256
#define MAX_RECORDS 100

// Struct to hold each frequency record
typedef struct 
{
    long long start_frequency;
    long long end_frequency;
    char service[MAX_SERVICE_LENGTH];
    char notes[MAX_NOTES_LENGTH];
} FrequencyRecord;

FrequencyRecord* read_table();
int lookup_record(FrequencyRecord*, long long);
long long get_long_from_float(float);

#endif