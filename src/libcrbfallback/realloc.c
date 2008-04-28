 #undef realloc

#include <config.h>
#include <sys/types.h>

void *realloc ();

void *
rpl_realloc (void* ptr, size_t n)
{
   if (n == 0)
   {
      n = 1;
      
      free (ptr);
      ptr = NULL;
   }

   return (ptr == NULL ? malloc (n) : realloc (ptr, n));
}
