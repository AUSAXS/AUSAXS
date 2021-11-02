SUBDIRS = lib tools test include doc

dist_noinst_DATA = COPYING.txt

dist-test: dist
	tar -z -x -f dsr-pdb*.tar.gz -C /tmp
	cd /tmp/dsr-pdb*; ./configure; make
