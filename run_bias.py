import MDAnalysis as mda
from MDAnalysis.analysis import pca, align, rms
import numpy as np
import os


def read_pdb(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None
    
fileout = open('cv_loss.csv', 'w')
fileout.write('time, cv_loss,\n')
start_pos = 'start_45_10'
start_rmsd = 9.123421959607505

ref_target = mda.Universe('target/human_AAC_C_state_1okc.pdb')
ref_non_target = mda.Universe('target/human_AAC_M_state_6gci.pdb')
i_run = 46

while start_rmsd >= 2:
    print(start_rmsd)
    n_run=1
    os.system(f'gmx grompp -f step7_production.mdp -o run_{i_run}_{n_run}.tpr -c {start_pos}.gro -p topol.top -n index.ndx -maxwarn 5')
    os.system(f'gmx mdrun -v -deffnm run_{i_run}_{n_run}')
    
    ref_run = mda.Universe(f'run_{i_run}_{n_run}.gro', f'run_{i_run}_{n_run}.xtc')
    
    R = mda.analysis.rms.RMSD(ref_run, ref_target, select='protein and name CA')
    R.run()
    
    current_rmsd = R.rmsd.T

    R = mda.analysis.rms.RMSD(ref_run, ref_non_target, select='protein and name CA')
    R.run()
    current_non_target_rmsd = 7-R.rmsd.T[2]
    
    cv_loss = current_rmsd[2]+current_non_target_rmsd
    
    min_current_rmsd = np.min(cv_loss)
    
    if min_current_rmsd < start_rmsd:
        
        start_rmsd = min_current_rmsd
        os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o start_{i_run}_{n_run}.gro -b {int(current_rmsd[:,np.argmin(cv_loss)][1])-1} -e {int(current_rmsd[:,np.argmin(cv_loss)][1])} < answers_out.txt')
        os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o result/run_{i_run}_{n_run}.xtc -b 0 -e {int(current_rmsd[:,np.argmin(cv_loss)][1])} < answers_out.txt')
        start_pos = f'start_{i_run}_{n_run}'
        fileout.write(f'{np.argmin(cv_loss)*2}, {start_rmsd},\n')
        i_run+=1
    else:
        
        while min_current_rmsd>=start_rmsd:
            
            print(start_rmsd)
            if n_run>20:
                os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o start_{i_run}_{n_run}.gro -b 249 -e 250 < answers_out.txt')
                os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o result/run_{i_run}_{n_run}.xtc -b 0 -e 250 < answers_out.txt')
                start_rmsd = min_current_rmsd
                start_pos = f'start_{i_run}_{n_run}'
                fileout.write(f'{250}, {start_rmsd},\n')
                i_run+=1
                break
            n_run+=1
            os.system(f'gmx grompp -f step7_production.mdp -o run_{i_run}_{n_run}.tpr -c {start_pos}.gro -p topol.top -n index.ndx -maxwarn 5')
            os.system(f'gmx mdrun -v -deffnm run_{i_run}_{n_run}')

            ref_run = mda.Universe(f'run_{i_run}_{n_run}.gro', f'run_{i_run}_{n_run}.xtc')

            R = mda.analysis.rms.RMSD(ref_run, ref_target, select='protein and name CA')
            R.run()

            current_rmsd = R.rmsd.T

            R = mda.analysis.rms.RMSD(ref_run, ref_non_target, select='protein and name CA')
            R.run()
            current_non_target_rmsd = 7-R.rmsd.T[2]

            cv_loss = current_rmsd[2]+current_non_target_rmsd

            min_current_rmsd = np.min(cv_loss)
            
            if min_current_rmsd < start_rmsd:

                start_rmsd = min_current_rmsd
                os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o start_{i_run}_{n_run}.gro -b {int(current_rmsd[:,np.argmin(cv_loss)][1])-1} -e {int(current_rmsd[:,np.argmin(cv_loss)][1])} < answers_out.txt')
                os.system(f'gmx trjconv -f run_{i_run}_{n_run}.xtc -s run_{i_run}_{n_run}.gro -o result/run_{i_run}_{n_run}.xtc -b 0 -e {int(current_rmsd[:,np.argmin(cv_loss)][1])} < answers_out.txt')
                
                start_pos = f'start_{i_run}_{n_run}'
                fileout.write(f'{np.argmin(cv_loss)*2}, {start_rmsd},\n')
                i_run+=1
                break
fileout.close()            