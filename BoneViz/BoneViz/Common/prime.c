#include "prime.h"
#include "set.h"

#define MAX_NUMBER 65536
#define MAX_NUM_PRIMES 65536

static int erato_table[MAX_NUM_PRIMES];
static int total_primes = 0;

static void 
fill_table()
{
  int table[MAX_NUMBER];
  int i;
  
  for (i=0; i<MAX_NUMBER; i++) table[i] = 0;

  for (i=2; i<MAX_NUMBER; i++) 
    if (table[i] == 0) {
      int j;
      
      for (j=2*i; j<MAX_NUMBER; j+=i) table[j] = 1;
    }
  
  for (i=2; i<MAX_NUMBER; i++)
    if (table[i] == 0) erato_table[total_primes++] = i;
}


boolean_t_new 
is_prime(int n)
{
  int i;

  if (total_primes == 0) fill_table();

  for (i=0; i<total_primes; i++) {
    if (n % erato_table[i] == 0) return FALSE;
    if (erato_table[i] * erato_table[i] > n) return TRUE;
  }

  return TRUE;
}

int 
next_prime(int n)
{
	int trueStuff = 1;	///to elim stupid compiler warning
	int result = 0;

	while(trueStuff){
		if (is_prime(n)){
			result = n;
			break;
		}
		else n++;
	}

	return n;
}
