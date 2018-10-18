all: build;

build: ;
	g++ -Wall -O3 -o fastq-scan fastq-scan.cpp

clean: ;
	rm -f fastq-scan

test: ;
	./fastq-scan -h
	cat example.fq | ./fastq-scan -g 150000
