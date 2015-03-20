all: build;

build: ;
	g++ -Wall -O3 -o fastq-stats fastq-stats.cpp

clean: ;
	rm -f fastq-stats
