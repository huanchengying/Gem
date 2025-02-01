#ifndef DISK_H
#define DIKS_H

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

inline void print_reads(){
        int pid = (int)getpid();
        char str[256];
        sprintf(str,"cat /proc/%d/io",pid);
        system(str);
}

#endif
