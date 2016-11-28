#include "../include/system.h"

void *poolptr;    // memory pool pointer
void *tmppool;    // temporary space pointer
int poolsize;     // size of memory pool
int tmpsize;      // size of temporary space
int curpool;      // current occupied space in pool
int curtmp;       // current occupied space in temp. space

// prints & flushes
void 
mprintf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
}

// prints error message and exits
void 
errexit(const char *format, ...) {
  va_list args;
  va_start(args, format);
  printf("***ERROR************\n");
  vprintf(format, args);
  printf("********************\n");
  va_end(args);
  exit(0);
}

// initializes memory pool
void
initMemoryPool(int ps, int ts){
  poolsize = ps;
  tmpsize  = ts;
  curpool  = 0;
  curtmp   = 0;
  if ( (poolptr  = malloc(ps))== NULL )
    errexit("Cannot allocate initial memory pool of %d bytes\n", ps);
  if ( (tmppool  = malloc(ts))== NULL )
    errexit("Cannot allocate initial memory pool of %d bytes\n", ts);
}


// alocates size bytes in memory pool
void *
smalloc(int size){
  void *retptr;
  if ( (size < 0) || (size > poolsize) )
    errexit("Reuqested block size (%d) out of allocation range ([0, %d])\n",
	    size,poolsize);
  if ( (curpool+size) > poolsize ){
     if ( (poolptr  = malloc(poolsize))== NULL )
       errexit("Cannot allocate additional memory pool\n");
     curpool = 0;
  }
  retptr   = poolptr;
  poolptr += size;
  curpool += size;
  return retptr;
}

// allocates size bytes of temporary space
void *
tmalloc(int size){
  void *retptr = tmppool;
  if ( (size < 0) || (size > tmpsize) )
    errexit("Reqested temporary block size (%d) out of allocation range ([0, %d])\n",
	    size,tmpsize);
  if ( (curtmp+size) > tmpsize )
    errexit("Temporary storage size (%d bytes) exceeded by  %d bytes (%d bytes requested)\n", tmpsize, 
	    curtmp+size-tmpsize, size);
  tmppool += size;
  curtmp  += size;
  return retptr;
}

// deallocates size bytes of temporary space which was allocated last
// and starts from ptr
void
tfree(void *ptr, int size){
  tmppool -= size;
  if ( tmppool != ptr )
    errexit("Somebody is trying to deallocate space (%d bytes) it doesn't own %d--%d\n",size, tmppool, ptr);
  curtmp -= size;
}

// deallocates size bytes of pool which was allocated last
// and starts from ptr
void
sfree(void *ptr, int size){
  poolptr -= size;
  if ( poolptr != ptr )
    errexit("Somebody is trying to deallocate space (%d bytes) it doesn't own %d--%d\n",size, poolptr, ptr);
  curpool -= size;
}

unsigned int
pow2(int x){
  unsigned int y = 1;
  return (y << x);
}

void printArray(int *list, int n, char *name){
  int i;
  mprintf("%s: ",name);
  for (i=0; i<n; i++)
    mprintf("%d ",list[i]);
  mprintf("\n");
}
