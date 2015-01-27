cd test_data
echo "../jobmanager.py -f ../iqtree_test_cmds.txt -c 16" | qsub -V -S /bin/bash -cwd -j y -r y -N iqtree_system_test -l zuseX -l cluster -pe threads 16 -q q.norm@zuse02  
