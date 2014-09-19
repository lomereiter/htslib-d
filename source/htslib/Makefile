all:
	rdmd --build-only -I.. -I../biod -L-L../c/htslib -L-lhts cram.d

optimized:
	rdmd --build-only --force --compiler=ldmd2 -I.. -O -release -inline -I../BioD -L-L../c/htslib -L-lhts cram.d
