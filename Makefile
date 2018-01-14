demultiplex_fastq:
	$(MAKE) -C src $@
	mkdir -p bin
	cp src/demultiplex_fastq bin

%:
	$(MAKE) -C src $@

