
#include "BuddyAllocator.h"

BuddyAllocator ba;


extern "C" void gatorInit( size_t bytes ) {
  auto myalloc = [] (size_t bytes) { void* ptr; cudaMallocManaged(&ptr,bytes); return ptr; };
  auto myfree  = [] (void* ptr) {cudaFree(ptr);};
  ba = BuddyAllocator( bytes , 1024 , myalloc , myfree );
}


extern "C" void gatorFinalize( ) {
  ba = BuddyAllocator();
}


extern "C" void* gatorAllocate( size_t bytes ) {
  return ba.allocate( bytes );
}


extern "C" void gatorDeallocate( void *ptr ) {
  ba.free( ptr );
}


