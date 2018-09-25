#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

void∗ Tokenize(void∗ rank) {
	long my rank = (long) rank;
	int count;
	int next = (my rank + 1) % thread count; char ∗fg rv;
	char my line[MAX]; char ∗my string;
	sem wait(&sems[my rank]);
	fg rv = fgets(my line, MAX, stdin); sem post(&sems[next]);
	while (fg rv != NULL) {
		printf("Thread %ld > my line = %s", my rank, my line);
		count = 0;
		my string = strtok(my line, " \t\n"); while ( my string != NULL ) {
		count++;
		printf("Thread %ld > string %d = %s\n", my rank, count,
		my string);
		my string = strtok(NULL, " \t\n");
	}
	sem wait(&sems[my rank]);
	fg rv = fgets(my line, MAX, stdin); sem post(&sems[next]);
	}
	return NULL;
} 