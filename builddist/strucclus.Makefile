CC     = cc -O3
OFILES = strucclus.o clusterlib.o 
LFILES = bioplib/fgetsany.o bioplib/array2.o bioplib/upstrncmp.o bioplib/GetWord.o \
         bioplib/OpenStdFiles.o bioplib/ReadPDB.o bioplib/ParseRes.c \
         bioplib/fsscanf.o bioplib/FreeStringList.o bioplib/StoreString.o \
         bioplib/FindNextResidue.o bioplib/WritePDB.o bioplib/chindex.o \
         bioplib/padterm.o bioplib/FindResidue.o bioplib/hash.o bioplib/PDBHeaderInfo.o \
         bioplib/BuildConect.o bioplib/throne.o bioplib/GetPDBChainLabels.o \
         bioplib/stringutil.o bioplib/prime.o bioplib/strcatalloc.o bioplib/PDB2Seq.o \
         bioplib/stringcat.o bioplib/IndexPDB.o bioplib/FindResidueSpec.o

all : strucclus


strucclus :  $(OFILES) $(LFILES)
	$(CC) -o strucclus $(OFILES) $(LFILES) -lm

.c.o  :
	$(CC) -o $@ -c $<

clean :
	\rm -f $(OFILES) $(LFILES)
