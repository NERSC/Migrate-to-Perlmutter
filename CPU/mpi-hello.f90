
program helloWorld
 implicit none
 include "mpif.h"
 integer :: myPE, numProcs, ierr
 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD, myPE, ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, ierr)
 print *, "Hello from Processor ", myPE 
 call MPI_FINALIZE(ierr)
end program helloWorld
