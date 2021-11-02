EXTRA_DIST = check_pdb_align_0.pdb        check_pdb_split.0.pdb\
check_pdb_align_1.pdb        check_pdb_split.1.pdb\
check_pdb_align_2.pdb        check_protein.pdb\
check_pdb_align_aligned.pdb  check_pdb_distance_output_0\
check_pdb_cat.pdb            check_pdb_distance_output_1\
check_pdb.pdb		     check_bonds.pdb


AM_CPPFLAGS = -I../include 
AM_CXXFLAGS = 

noinst_PROGRAMS = check_protein check_pdb check_transform check_refine_alignment check_bonds

check_protein_SOURCES = check_protein.cc
check_protein_LDADD = ../lib/libdsrpdb.la 
check_pdb_SOURCES = check_pdb.cc
check_pdb_LDADD = ../lib/libdsrpdb.la 
check_transform_SOURCES = check_transform.cc
check_transform_LDADD = ../lib/libdsrpdb.la 
check_refine_alignment_SOURCES = check_refine_alignment.cc
check_refine_alignment_LDADD = ../lib/libdsrpdb.la 
check_bonds_SOURCES = check_bonds.cc
check_bonds_LDADD = ../lib/libdsrpdb.la 

TEMPTAG=check_out
DIFF = diff -w -b -B

test: check test_cleanup test_transform test_reading test_cat test_split test_align  test_distance test_bonds test_cleanup

test_cleanup:
	@rm -f *$(TEMPTAG)*

test_transform:
	@echo -n Checking Transform...
	@./check_transform check_protein.pdb > /dev/null
	@echo done.

test_reading:
	@echo -n Checking basic reading and writing facilities...
	@./check_protein check_protein.pdb check_protein.$(TEMPTAG).pdb > /dev/null
	@$(DIFF) check_protein.pdb check_protein.$(TEMPTAG).pdb
	@./check_pdb check_pdb.pdb check_pdb.$(TEMPTAG).pdb > /dev/null
	@$(DIFF) check_pdb.pdb check_pdb.$(TEMPTAG).pdb
	@echo done.

test_cat:
	@echo -n Checking pdb_cat...
	@../tools/pdb_cat check_protein.pdb check_protein.pdb check_pdb_cat.$(TEMPTAG).pdb
	@$(DIFF) check_pdb_cat.pdb check_pdb_cat.$(TEMPTAG).pdb
	@echo done.

test_split:
	@echo -n Checking pdb_split...
	@../tools/pdb_split check_pdb_cat.pdb check_pdb_split.%d.pdb > /dev/null
	@../tools/pdb_cat  check_pdb_split.*.pdb check_pdb_cat.$(TEMPTAG).pdb
	@$(DIFF) check_pdb_cat.pdb check_pdb_cat.$(TEMPTAG).pdb
	@echo done.

test_align:
	@echo -n Checking pdb_align...
	@../tools/pdb_align check_pdb_align_1.pdb check_pdb_align_0.pdb -a check_pdb_align_aligned.$(TEMPTAG).pdb > /dev/null
	@$(DIFF) check_pdb_align_aligned.pdb  check_pdb_align_aligned.$(TEMPTAG).pdb 
	@echo done.

test_distance:
	@echo -n Checking pdb_distance...
	@../tools/pdb_distance check_pdb_align_0.pdb check_pdb_align_0.pdb -d -c > check_pdb_distance_0.$(TEMPTAG)
	@$(DIFF) check_pdb_distance_output_0 check_pdb_distance_0.$(TEMPTAG)
	@../tools/pdb_distance check_pdb_align_0.pdb check_pdb_align_2.pdb -d -c > check_pdb_distance_1.$(TEMPTAG)
	@$(DIFF) check_pdb_distance_output_1 check_pdb_distance_1.$(TEMPTAG)
	@echo done.

test_bonds:
	@echo -n Checking bond interators...
	@./check_bonds
	@echo done.
