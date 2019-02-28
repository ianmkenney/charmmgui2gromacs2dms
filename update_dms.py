import sqlite3
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mda

def update_positions_and_velocities(gro, dmsconnection, mapping):
    """Iterates over the gro atoms and updates the positions and velocities of the 
    atoms in the DMS file.

    Args:
        gro (Universe): Universe from the gro file
        dmsconnection (sqlite connection): Connection to the DMS file
        mapping (dict): map of the indices from the gro to the DMS format

    Returns:
        None
    """
    
    c = dmsconnection.cursor()

    for atom in gro.atoms:
        dms_index = mapping[str(atom.index)]
        x, y, z = atom.position
        vx, vy, vz = atom.velocity
        
        c.execute('''select * from particle where id == {0}'''.format(dms_index))
        Q = c.fetchone()

        if Q[2] == atom.name:
            c.execute('''UPDATE PARTICLE 
                         SET x = {x}, 
                             y = {y}, 
                             z = {z}, 
                             vx = {vx}, 
                             vy = {vy}, 
                             vz = {vz} 
                        WHERE id == {index}'''.format(x=x, y=y, z=z, 
                                                      vx=vx, vy=vy, vz=vz, 
                                                      index=dms_index))
        else:
            raise ValueError

def update_pbc_box(gro, dmsconnection):
    """Update the PBC box in the DMS from the gro file.

    Args:
        gro (Universe): Universe from the gro file
        dmsconnection (sqlite connection): Connection to the DMS file

    Returns:
        None
    """
    
    c = dmsconnection.cursor()
    dims = list("xyz")
    for i, d in enumerate(dims):
        dmsid = i+1
        c.execute("""UPDATE global_cell
                     SET {dim} = {dimvalue}
                     WHERE id == {dmsid}""".format(dim=d, dimvalue=gro.dimensions[i], dmsid=dmsid))

if __name__ == "__main__":

    import sys
    import json


    if "-h" in sys.argv:
        help_msg = "example usage: update_dms.py map.json pos_and_vel.gro mutate.dms"
        print(help_msg)
        exit(1)

    with open(sys.argv[1],'r') as f:
        mapping = json.load(f)

    gro = mda.Universe(sys.argv[2])
    dmsconnection = sqlite3.connect(sys.argv[3])

    print("Updating DMS positions and velocities")
    update_positions_and_velocities(gro, dmsconnection, mapping)
    print("Updating DMS box dimensions")
    update_pbc_box(gro, dmsconnection)

    print("Committing changes")
    dmsconnection.commit()
    dmsconnection.close()
