OPTS   = -O3 -L$(HOME)/lib -I$(HOME)/include
OFILES = cluster.o clusterlib.o
IFILES = cluster.h
TARGET = cluster
LIBS   = -lgen -lm


$(TARGET) : $(OFILES)
	$(CC) $(OPTS) -o $@ $(OFILES) $(LIBS)

.c.o : $(IFILES)
	$(CC) $(OPTS) -c -o $@ $<

clean :
	\rm -f $(OFILES)

