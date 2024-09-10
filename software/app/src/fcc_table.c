#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fcc_table.h"


void read_csv(const char *filename, FrequencyRecord *records, int *count) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    int i = 0;

    // Skip the header line
    fgets(line, sizeof(line), file);

    while (fgets(line, sizeof(line), file) && i < MAX_RECORDS) {
        // Parse the CSV line
        char *token;

        // Start frequency
        token = strtok(line, ",");
        records[i].start_frequency = atoll(token);

        // End frequency
        token = strtok(NULL, ",");
        records[i].end_frequency = atoll(token);

        // Service/Use
        token = strtok(NULL, ",");
        strncpy(records[i].service, token, MAX_SERVICE_LENGTH);
        records[i].service[MAX_SERVICE_LENGTH - 1] = '\0';  // Ensure null-termination

        // Notes
        token = strtok(NULL, ",");
        strncpy(records[i].notes, token, MAX_NOTES_LENGTH);
        records[i].notes[MAX_NOTES_LENGTH - 1] = '\0';  // Ensure null-termination

        i++;
    }

    *count = i;

    fclose(file);
}

long long get_long_from_float(float frequency_float)
{
    return (long)frequency_float;      
}

int lookup_record(FrequencyRecord* precord_table, long long freq_to_lookup)
{
    int index = 0;
    for (index = 0; index < MAX_RECORDS; index++)
    {
        if (precord_table[index].start_frequency < freq_to_lookup && 
            precord_table[index].end_frequency > freq_to_lookup)
            return index;
    }

    return -1;
}

FrequencyRecord* read_table() {
    //FrequencyRecord records[MAX_RECORDS];
    FrequencyRecord* records = malloc(MAX_RECORDS * sizeof(FrequencyRecord));
    int count = 0;

    // Read the CSV file into the array of structs
    read_csv("fcc_table.csv", records, &count);
    
    return records;
}
