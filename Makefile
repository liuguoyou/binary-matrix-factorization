#TEST
all:
	make -C src
	make -C doc

clean:
	make -C src clean
	make -C doc clean
