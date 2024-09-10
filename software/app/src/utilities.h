#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <termios.h>
#include <fcntl.h>
#include <hackrf.h>

int kbhit(void);
int parse_int(char *s, int *const value);
int parse_u32(char *s, uint32_t *const value);
int parse_u32_range(char *s, uint32_t *const value_min, uint32_t *const value_max);

#endif