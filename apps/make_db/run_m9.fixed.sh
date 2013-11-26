export LD_LIBRARY_PATH=~/perm-je-0.9.6/lib
#./make_db_table.32bit.20mer -l -i m9.f.lst -o /local/ramfs/m9.0904.16.db -k 20 -s 650 
./make_db_table.16bit.20mer -l -i m9.f.lst -o /local/ramfs/m9.20mer.16bit.g1000.db -k 20 -s 450 -g 1000 -m /usr/mic/post1/metagenomics/runtime_inputs/09042013/numeric_ranks -f /usr/mic/post1/metagenomics/runtime_inputs/09042013/m9.32To16.map

