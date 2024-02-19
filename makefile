CC = gcc
CFLAGS = -I. -I/opt/homebrew/opt/gsl/include
LDFLAGS = -L/opt/homebrew/opt/gsl/lib -lgsl -lgslcblas -lm

BIN = cluster computeMSE computeNumSzReads createFastq generateStats NtoMinQV computeOptimalRho computeSvd compressQval decompressQval
all: $(BIN)

cluster: cluster.c
	$(CC) $(CFLAGS) -o cluster cluster.c $(LDFLAGS)

computeMSE: computeMSE.c
	$(CC) $(CFLAGS) -o computeMSE computeMSE.c $(LDFLAGS)

computeNumSzReads: computeNumSzReads.c
	$(CC) $(CFLAGS) -o computeNumSzReads computeNumSzReads.c

createFastq: createFastq.c
	$(CC) $(CFLAGS) -o createFastq createFastq.c

generateStats: generateStats.c
	$(CC) $(CFLAGS) -o generateStats generateStats.c

NtoMinQV: NtoMinQV.c
	$(CC) $(CFLAGS) -o NtoMinQV NtoMinQV.c

computeOptimalRho: computeOptimalRho.c
	$(CC) $(CFLAGS) -o computeOptimalRho computeOptimalRho.c $(LDFLAGS)

computeSvd: computeSvd.c
	$(CC) $(CFLAGS) -o computeSvd computeSvd.c $(LDFLAGS)

compressQval: compressQval.c
	$(CC) $(CFLAGS) -o compressQval compressQval.c mtwist.c randistrs.c $(LDFLAGS)

decompressQval: decompressQval.c
	$(CC) $(CFLAGS) -o decompressQval decompressQval.c mtwist.c randistrs.c $(LDFLAGS)

clean:
	$(RM) $(BIN) *.o
