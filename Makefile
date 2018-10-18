all: build;

build: ;
	g++ -Wall -O3 -o fastq-scan fastq-scan.cpp

clean: ;
	rm -f fastq-scan
