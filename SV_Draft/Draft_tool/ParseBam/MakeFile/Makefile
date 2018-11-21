SharedLibraryPath=../../../../ShareLibrary/



BAMTOOLS=../../../../bamtools/install/include/bamtools
BAMTOOLS_LD=../../../../bamtools/install/lib

CC = g++
INC_DIR += $(SharedLibraryPath) #Identify include path
INC_DIR += $(BAMTOOLS) #Identify include path
CFLAGS = -O3 -Wall -std=c++0x -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
LIBS += -lz

OBJ = ./obj/

all: ParseBam

clsbasealgorithm.o : $(SharedLibraryPath)clsbasealgorithm.cpp $(SharedLibraryPath)clsbasealgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbasealgorithm.cpp

clsfastareader.o: $(SharedLibraryPath)clsfastareader.cpp $(SharedLibraryPath)clsfastareader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastareader.cpp

clsfastqreader.o: $(SharedLibraryPath)clsfastqreader.cpp $(SharedLibraryPath)clsfastqreader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastqreader.cpp


clsbwa.o : $(SharedLibraryPath)clsbwa.cpp $(SharedLibraryPath)clsbwa.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbwa.cpp

clsvcf1000genome.o : $(SharedLibraryPath)clsvcf1000genome.cpp $(SharedLibraryPath)clsvcf1000genome.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsvcf1000genome.cpp

clsreadconfigini.o : $(SharedLibraryPath)clsreadconfigini.cpp $(SharedLibraryPath)clsreadconfigini.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsreadconfigini.cpp

clsmuscle.o : $(SharedLibraryPath)clsmuscle.cpp $(SharedLibraryPath)clsmuscle.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsmuscle.cpp

clsvelvet.o : $(SharedLibraryPath)clsvelvet.cpp $(SharedLibraryPath)clsvelvet.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsvelvet.cpp

clsblast.o : $(SharedLibraryPath)clsblast.cpp $(SharedLibraryPath)clsblast.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsblast.cpp

clskmeralgorithm.o : $(SharedLibraryPath)clskmeralgorithm.cpp $(SharedLibraryPath)clskmeralgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clskmeralgorithm.cpp

clsparsebam.o: ../ParseBam/clsparsebam.cpp ../ParseBam/clsparsebam.h
	$(CC) $(CFLAGS) -c ../ParseBam/clsparsebam.cpp

clscomparison.o: ../ParseBam/clscomparison.cpp ../ParseBam/clscomparison.h
	$(CC) $(CFLAGS) -c ../ParseBam/clscomparison.cpp

clsdebug.o: ../ParseBam/clsdebug.cpp ../ParseBam/clsdebug.h
	$(CC) $(CFLAGS) -c ../ParseBam/clsdebug.cpp

clsdrawimage.o: ../ParseBam/clsdrawimage.cpp ../ParseBam/clsdrawimage.h
	$(CC) $(CFLAGS) -c ../ParseBam/clsdrawimage.cpp

clsconfig.o: ../ParseBam/clsconfig.cpp ../ParseBam/clsconfig.h
	$(CC) $(CFLAGS) -c ../ParseBam/clsconfig.cpp

main.o: ../ParseBam/main.cpp 
	$(CC) $(CFLAGS) -c ../ParseBam/main.cpp

ParseBam: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsvcf1000genome.o clsreadconfigini.o clsmuscle.o clsvelvet.o clsblast.o clskmeralgorithm.o clsparsebam.o clscomparison.o clsdebug.o clsdrawimage.o clsconfig.o main.o 
	$(CC) $(CFLAGS) -o ParseBam clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsvcf1000genome.o clsreadconfigini.o clsmuscle.o clsvelvet.o clsblast.o clskmeralgorithm.o clsparsebam.o clscomparison.o clsdebug.o clsdrawimage.o clsconfig.o main.o -lz -lbamtools


RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


