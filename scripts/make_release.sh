#!/bin/sh
#  set version number
verno=1.2.2a

reldir=../../LMAT-$verno

#create directories

if [ ! -d "$reldir" ] 
then

    mkdir $reldir
    mkdir $reldir/src
    mkdir $reldir/src/kmerdb
    mkdir $reldir/src/kmerdb/lib
    mkdir $reldir/include
    mkdir $reldir/doc
    mkdir $reldir/lib
    mkdir $reldir/bin
    mkdir $reldir/runtime_inputs
    mkdir $reldir/example
    mkdir $reldir/third-party
fi

#copy doc, set version number

cp ../lmat-doc.tex $reldir/doc

sed -i "s/VERSIONNO/$verno/g" $reldir/doc/lmat-doc.tex

pushd $reldir/doc

latex lmat-doc.tex 
latex lmat-doc.tex 
latex2html -split 0 lmat-doc.tex 

rm lmat-doc/WARNINGS

links -dump -no-numbering -no-references lmat-doc/index.html > lmat-doc.txt

#cleanup

for ext in tex aux log out toc dvi
do
    rm $reldir/doc/lmat-doc.$ext
done

popd

#copy files

cp release-templates/README.src $reldir/README


cp release-templates/Makefile.inc.src $reldir/Makefile.inc
cp release-templates/Makefile.src $reldir/Makefile
cp release-templates/Makefile.src.kmerdb $reldir/src/kmerdb/Makefile
cp release-templates/Makefile.src.src $reldir/src/Makefile

cp release-templates/third-party/* $reldir/third-party

# if we need to refer to the version number...
#sed -i "s/VERSIONNO/$verno/g" $reldir/Makefile.inc

# copy source files

#edit lists here
appsrc=( apps/read_label preprocessing/kmerPrefixCounter apps/make_db/make_db_table preprocessing/tax_histo_new_fmt apps/rand_read_label apps/content_summ apps/gene_label rmodel/countTaxidFrequency)

#for now this list needs to be manually sync'ed with release-templates/Makefile.src.src

## these should match up with contents of wrapper scripts
appdest=( read_label kmerPrefixCounter make_db_table tax_histo rand_read_label content_summ gene_label frequency_counter )

for i in `seq 0 7`
do
echo "copy: ${appsrc[$i]}"
cp ../${appsrc[$i]}.cpp $reldir/src/${appdest[$i]}.cpp
echo "${appdest[$i]}:" '$(METAG_LIB)' ${appdest[$i]}.cpp >> $reldir/src/Makefile
echo -e "\t"'$(CXX) -std=c++0x $(CXXFLAGS) $(LDFLAGS) $@.cpp $(LIBS) -o $@' >> $reldir/src/Makefile
echo -e "\t""cp ${appdest[$i]} ../bin" >> $reldir/src/Makefile
echo -e "\n"  >> $reldir/src/Makefile 
done


# copy kmerdb srcs
#pwd
for n in Utils.hpp KmerNode.hpp KmerFileMetaData.hpp TaxTree.hpp TaxNode.hpp TaxTable.hpp all_headers.hpp TaxNodeStat.hpp SortedDb.hpp Encoder.hpp KmerIterator.hpp StopWatch.hpp all_headers.hpp metag_typedefs.hpp
do
    echo "copy $n"
    cp ../src/kmerdb/$n $reldir/src/kmerdb
done


# assumes cpp extension
for n in Utils KmerFileMetaData TaxTree TaxNode TaxTable SortedDb
do
#    echo $n
    cp ../src/kmerdb/$n.cpp $reldir/src/kmerdb
done

# copy include

## I don't think kencode.hpp is used anymore
for n in metag_const.h perm.h tid_checks.hpp
do
cp ../include/$n $reldir/include
done



#copy preprocessing scripts
for n in merge_fastq_reads_with_N_separator.pl Tid16_getMapping.py build_tid_numeric_rank_table.py build_species_level_map.py Tid16_get32BitTaxIDs.py get_gi_numbers.py build_header_table.py
do
cp ../preprocessing/$n $reldir/bin
done

#copy nclass scripts

for n in run_lmat.sh tolineage.py combine_fast.pl losummary_fast.pl losummary_fast_mc.sh pull_reads_mc.sh pull_reads.pl build_taxid_lst.sh build_taxid_lst.pl
do
    if [ ! -f ../nclass/$n ] ; then
	echo $n not found
    else  
	cp ../nclass/$n $reldir/bin
    fi
done 

#copy rmodel scripts

for n in gen_rand_mod.sh merge_cnts.py
do
cp ../rmodel/$n  $reldir/bin
done

#copy download script

cp get_db.sh $reldir/bin

#copy example file from the server

cp /usr/mic/post1/metagenomics/metagenomes/synthetic/example.tgz $reldir/example

# copy runtime inputs (?)

cd ../..

tar -czvf LMAT-$verno.tgz LMAT-$verno/
mv LMAT-$verno.tgz /usr/mic/post1/metagenomics/tarballs

exit 0

readmefile=$reldir/runtime_inputs/README
echo "DESCRIPTION OF RUNTIME INPUT FILES" > $readmefile
echo "--------------------------------------" >> $readmefile
echo "ncbi_taxonomy.segment.dat - NCBI taxonomy tree in LMAT readable format" >> $readmefile
cp ../runtime_inputs/ncbi_taxonomy.segment.dat $reldir/runtime_inputs/
echo "ncbi_taxonomy.segment.dat.nohl - NCBI taxonomy tree in LMAT readable format, cuts out most of the human lineage to speed up processing" >> $readmefile
cp ../runtime_inputs/ncbi_taxonomy.segment.dat.nohl $reldir/runtime_inputs/
echo "ncbi_taxonomy.segment.dat - Human readable NCBI taxonomy lineage" >> $readmefile
cp ../runtime_inputs/ncbi_taxonomy_rank.segment.txt $reldir/runtime_inputs/
echo "ncbi_taxid_to_rank.txt - Report NCBI taxonomy rank for each NCBI taxonomy node" >> $readmefile
cp ../runtime_inputs/ncbi_taxid_to_rank.txt $reldir/runtime_inputs/
echo "depth_for_ncbi_taxonomy.segment.dat - Report NCBI taxonomy tree depth for each NCBI taxonomy node" >> $readmefile
cp ../runtime_inputs/depth_for_ncbi_taxonomy.segment.dat $reldir/runtime_inputs/
echo "m9_32To16BitMap.txt - Mapping of NCBI taxonomy IDs to 16-bit encoding (for compression)" >> $readmefile
cp ../runtime_inputs/m9_32To16BitMap.txt $reldir/runtime_inputs/
echo "supsect_m9_genomes.txt - User identified genomes in the database with possible contaminants, require higher threshold in final content summarization step" >> $readmefile
cp ../runtime_inputs/supsect_m9_genomes.txt $reldir/runtime_inputs/
echo "null.bin.10.m9.db.500.*.4800000.rl_output.rand_lst.gz - Random null models" >> $readmefile
cp ../runtime_inputs/null.bin.10.m9.db.500.*.4800000.rl_output.rand_lst.gz $reldir/runtime_inputs/
echo "null_lst.txt - List of null models to be used by LMAT and their associated k-mer count (read length - k - 1)" >> $readmefile
cp ../runtime_inputs/null_lst.txt $reldir/runtime_inputs/
echo "tcnt.m9.20.tax_histo - counts the number of distinct k-mers associated with each NCBI taxonomy ID" >> $readmefile
cp ../runtime_inputs/tcnt.m9.20.tax_histo $reldir/runtime_inputs/
echo "gn_ref2.txt.gz - Gene database" >> $readmefile
cp ../runtime_inputs/gn_ref2.txt.gz $reldir/runtime_inputs/
echo "low_numid_plasmids.txt - Identifies the plasmids in NCBI that are not associated with an organism taxonomy ID" >> $readmefile
cp ../runtime_inputs/low_numid_plasmids.txt $reldir/runtime_inputs/
