#include "int32_sort.h"
#include "uint32_sort.h"

/* can save time by vectorizing xor loops */
/* can save time by integrating xor loops with int32_sort */

void uint32_sort(uint32_t *x,long long n)
{
  long long j;
   // printf("x antes\n");
   //    for (int i = 0; i < n; ++i)
   //     {
   //       printf("%d ", x[i]);
   //     }
   //     printf("\n");

  for (j = 0;j < n;++j) x[j] ^= 0x80000000;

  	 // printf("x antes\n");
    //   for (int i = 0; i < n; ++i)
    //    {
    //      printf("%d ", x[i]);
    //    }
    //    printf("\n");
       
  int32_sort((int32_t *) x,n);
  for (j = 0;j < n;++j) x[j] ^= 0x80000000;
}
