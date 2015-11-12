#
# This is an example Makefile for a countwords program.  This
# program uses both the scanner module and a counter module.
# Typing 'make' or 'make count' will create the executable file.
#

# define some Makefile variables for the compiler and compiler flags
# to use Makefile variables later in the Makefile: $()
#
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#
# for C++ define  CC = g++
CC = icpc
CFLAGS  = -Wall

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: incr_pagerank

# To create the executable file count we need the object files
# countwords.o, counter.o, and scanner.o:
#
incr_pagerank:  sparse_vector.o sparse_matrix_csc.o incr_pagerank.o 
	$(CC) $(CFLAGS) -o incr_pagerank sparse_vector.o incr_pagerank.o sparse_matrix_csc.o

	

# To create the object file countwords.o, we need the source
# files countwords.c, scanner.h, and counter.h:
#
sparse_vector.o:  sparse_vector.cc 
	$(CC) $(CFLAGS) -c sparse_vector.cc

# To create the object file counter.o, we need the source files
# counter.c and counter.h:
#
incr_pagerank.o:  incr_pagerank.cc 
	$(CC) $(CFLAGS) -c incr_pagerank.cc

# To create the object file scanner.o, we need the source files
# scanner.c and scanner.h:
#
sparse_matrix_csc.o:  sparse_matrix_csc.cc 
	$(CC) $(CFLAGS) -c sparse_matrix_csc.cc

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) incr_pagerank *.o *~