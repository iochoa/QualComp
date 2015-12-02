CC = gcc
CFLAGS = -I.
LDFLAGS = -lgsl -lgslcblas -lm

BIN = cluster computeMSE computeNumSzReads createFastq generateStats NtoMinQV computeOptimalRho computeSvd compressQval decompressQval
all: $(BIN)

cluster: cluster.c
	$(CC) -o cluster cluster.c $(LDFLAGS)

computeMSE: computeMSE.c
	$(CC) -o computeMSE computeMSE.c $(LDFLAGS)

computeNumSzReads: computeNumSzReads.c
	$(CC) -o computeNumSzReads computeNumSzReads.c

createFastq: createFastq.c
	$(CC) -o createFastq createFastq.c

generateStats: generateStats.c
	$(CC) -o generateStats generateStats.c

NtoMinQV: NtoMinQV.c
	$(CC) -o NtoMinQV NtoMinQV.c

computeOptimalRho: computeOptimalRho.c
	$(CC) -o computeOptimalRho computeOptimalRho.c $(LDFLAGS)

computeSvd: computeSvd.c
	$(CC) -o computeSvd computeSvd.c $(LDFLAGS)

compressQval: compressQval.c
	$(CC) -I ./  -o compressQval compressQval.c mtwist.c randistrs.c $(LDFLAGS)

decompressQval: decompressQval.c
	$(CC) -I ./ -o decompressQval  decompressQval.c mtwist.c randistrs.c $(LDFLAGS)

clean:
	$(RM) $(BIN) *.o
