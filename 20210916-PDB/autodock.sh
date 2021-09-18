pythonsh ../../Utilities24/prepare_receptor4.py -r ./5Z15.pdb -A hydrogens -U waters
pythonsh ../../Utilities24/prepare_ligand4.py -l ./DB01611.pdb -A hydrogens -U waters

# pythonsh ../Utilities24/prepare_dpf4.py -l DB01118.pdbqt -r 1GQ4.pdbqt
# pythonsh ../Utilities24/prepare_gpf4.py -l DB01118.pdbqt -r 1GQ4.pdbqt

# autogrid4 –p 1GQ4.gpf
# autodock4 –p DB01118_1GQ4.dpf
~/autodock_vina_1_1_2_linux_x86/bin/vina --receptor ./P07550/1GQ4.pdbqt --ligand DB01118.pdbqt --exhaustiveness 30 --center_x 0 --center_y 0 --center_z 0 --size_x 30 --size_y 30 --size_z 30 --log 1GQ4_DB01118_1GQ4.txt