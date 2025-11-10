import re,random, math,argparse
import numpy as np
from collections import defaultdict

def euler_rotate(xyz, axis="z"):
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    theta = random.randint(0, 360) * (2 * math.pi / 360)
    if axis == "x":
        x1 = x
        y1 = y * np.cos(theta) - z * np.sin(theta)
        z1 = y * np.sin(theta) + z * np.cos(theta)
    elif axis == "y":
        x1 = x * np.cos(theta) + z * np.sin(theta)
        y1 = y
        z1 = -x * np.sin(theta) + z * np.cos(theta)
    elif axis == "z":
        x1 = x * np.cos(theta) - y * np.sin(theta)
        y1 = x * np.sin(theta) + y * np.cos(theta)
        z1 = z
    xyz1 = np.zeros(np.shape(xyz))
    xyz1[:, 0] = x1
    xyz1[:, 1] = y1
    xyz1[:, 2] = z1
    return xyz1

def parse_lammps_data(content):
    sections = defaultdict(list)
    current_section = None
    for line in content.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("LAMMPS data file"):
            continue
        match = re.match(r'^(Masses|Pair Coeffs|Bond Coeffs|Angle Coeffs|Dihedral Coeffs|Improper Coeffs|Atoms|Bonds|Angles|Dihedrals|Impropers)$', stripped)
        if match:
            current_section = match.group()
            continue
        if current_section and re.match(r'^\d', stripped):
            sections[current_section].append(stripped)
    return sections

def renumber_lines(lines, id_offset, type_offset, atom_offset):
    new_lines = []
    for line in lines:
        parts = line.split()
        new_id = int(parts[0]) + id_offset
        new_type = int(parts[1]) + type_offset
        new_parts = [str(new_id), str(new_type)] + [str(int(p) + atom_offset) for p in parts[2:]]
        new_lines.append(" ".join(new_parts))
    return new_lines

def random_translate_atoms(atom_lines, dx, dy, dz):
    translated_lines = []
    for line in atom_lines:
        parts = line.split()
        x, y, z = map(float, parts[4:7])
        new_coords = [str(x + dx), str(y + dy), str(z + dz)]
        translated_line = " ".join(parts[:4] + new_coords + parts[7:])
        translated_lines.append(translated_line)
    return translated_lines

def generate_random_position(box_size, padding=3.0):
    return np.random.uniform(low=padding, high=box_size - padding, size=3)

def merge_multiple_lammps_files(lmp_files, mol_counts, output_path, box_size=60.0):
    assert len(lmp_files) == len(mol_counts), "Each LAMMPS file must have a corresponding molecule count."

    cumulative = defaultdict(list)
    total_counts = {'Atoms': 0, 'Bonds': 0, 'Angles': 0, 'Dihedrals': 0, 'Impropers': 0,
                    'atom types': 0, 'bond types': 0, 'angle types': 0, 'dihedral types': 0, 'improper types': 0}
    atom_offset = bond_offset = angle_offset = dihedral_offset = improper_offset = 0
    atom_type_offset = bond_type_offset = angle_type_offset = dihedral_type_offset = improper_type_offset = 0
    mol_id = 1

    for idx, (lmp_file, mol_count) in enumerate(zip(lmp_files, mol_counts)):
        with open(lmp_file, 'r') as f:
            content = f.read()

        parsed = parse_lammps_data(content)

        # Masses
        new_masses = [
            f"{int(line.split()[0]) + atom_type_offset} {line.split()[1]}"
            for line in parsed.get('Masses', [])
        ]
        cumulative['Masses'].extend(new_masses)

        # Coeffs
        for key, offset in zip(
            ['Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs'],
            [atom_type_offset, bond_type_offset, angle_type_offset, dihedral_type_offset, improper_type_offset]
        ):
            new_block = [
                f"{int(line.split()[0]) + offset} " + " ".join(line.split()[1:])
                for line in parsed.get(key, [])
            ]
            cumulative[key].extend(new_block)

        for i in range(mol_count):
            dx, dy, dz = generate_random_position(box_size)
            atom_lines = [
                f"{int(p.split()[0]) + atom_offset} {mol_id} {int(p.split()[2]) + atom_type_offset} " + " ".join(p.split()[3:])
                for p in parsed.get('Atoms', [])
            ]
            atom_lines = random_translate_atoms(atom_lines, dx, dy, dz)
            cumulative['Atoms'].extend(atom_lines)

            cumulative['Bonds'].extend(renumber_lines(parsed.get('Bonds', []), bond_offset, bond_type_offset, atom_offset))
            cumulative['Angles'].extend(renumber_lines(parsed.get('Angles', []), angle_offset, angle_type_offset, atom_offset))
            cumulative['Dihedrals'].extend(renumber_lines(parsed.get('Dihedrals', []), dihedral_offset, dihedral_type_offset, atom_offset))
            cumulative['Impropers'].extend(renumber_lines(parsed.get('Impropers', []), improper_offset, improper_type_offset, atom_offset))

            atom_offset += len(parsed.get('Atoms', []))
            bond_offset += len(parsed.get('Bonds', []))
            angle_offset += len(parsed.get('Angles', []))
            dihedral_offset += len(parsed.get('Dihedrals', []))
            improper_offset += len(parsed.get('Impropers', []))
            mol_id += 1

        total_counts['atom types'] += len(parsed.get('Masses', []))
        total_counts['bond types'] += len(parsed.get('Bond Coeffs', []))
        total_counts['angle types'] += len(parsed.get('Angle Coeffs', []))
        total_counts['dihedral types'] += len(parsed.get('Dihedral Coeffs', []))
        total_counts['improper types'] += len(parsed.get('Improper Coeffs', []))

        atom_type_offset = total_counts['atom types']
        bond_type_offset = total_counts['bond types']
        angle_type_offset = total_counts['angle types']
        dihedral_type_offset = total_counts['dihedral types']
        improper_type_offset = total_counts['improper types']

    with open(output_path, 'w') as f:
        f.write("LAMMPS data file merged from multiple sources\n\n")
        f.write(f"{len(cumulative['Atoms'])} atoms\n")
        f.write(f"{len(cumulative['Bonds'])} bonds\n")
        f.write(f"{len(cumulative['Angles'])} angles\n")
        f.write(f"{len(cumulative['Dihedrals'])} dihedrals\n")
        f.write(f"{len(cumulative['Impropers'])} impropers\n\n")

        f.write(f"{total_counts['atom types']} atom types\n")
        f.write(f"{total_counts['bond types']} bond types\n")
        f.write(f"{total_counts['angle types']} angle types\n")
        f.write(f"{total_counts['dihedral types']} dihedral types\n")
        f.write(f"{total_counts['improper types']} improper types\n\n")

        f.write(f"0.0 {box_size} xlo xhi\n0.0 {box_size} ylo yhi\n0.0 {box_size} zlo zhi\n\n")

        for section in ['Masses', 'Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Impropers']:
            if cumulative[section]:
                f.write(f"{section}\n\n")
                for line in cumulative[section]:
                    f.write(line + "\n")
                f.write("\n")

    with open('settings.lmp', 'w') as f:
        for section in ['Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs']:
            if cumulative[section]:
                label = '_'.join(section.split(' ')).lower()[:-1]
                for line in cumulative[section]:
                    if section == 'Pair Coeffs':
                        f.write(f'{label}' + '\t' + f'{line.split()[0]}' + '\t' + line + "\n")
                    else:
                        f.write(f'{label}' + '\t' + line + "\n")
                f.write("\n")

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='''This script reads a .lmp file output by OPLS and creates a box and writes an init file for it.                                                                                                                                          
    Input: .lmp file created by OPLS                                                                                                                                                                                                  
    Output: .init file and new .lmp/data file.''')

    # Required Arguments
    parser.add_argument('lmp', type=str, help='.lmp file for the desired molecule provided by ligpargen.')

    # Optional Arguments
    parser.add_argument('-N', dest='N', type=str, default='10', help='Number of molecules to be added to the simulation box. ')
    parser.add_argument('-O', dest='output', type=str, default='data.lmp', help='Name of .lmp file to be produced. ')
    parser.add_argument('-L', dest='length', type=int, default=30.0, help='Length of the box in Angstroms. Default is 30.0 Angstroms. ')
    
    # Parse Arguments
    args = parser.parse_args()

    lmp_files = args.lmp.split()
    mol_counts = [int(n) for n in args.N.split()]
    assert len(lmp_files) == len(mol_counts), "Number of .lmp files must match the number of molecule counts provided."
    
    merge_multiple_lammps_files(lmp_files, mol_counts, args.output, args.length)
