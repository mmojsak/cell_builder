print('--- importing packages ---')
print()

#########################################################################################
### IMPORT PACKAGES
# Atomic Sumilation Environment
import ase
from ase import Atoms
from ase.io import read, write
# Mathematics
import numpy as np
from numpy import arccos, radians
import skspatial
from skspatial.objects import Line, Point
# Object manipulation
import copy
# User interface
import sys
import warnings # used to suppress CASTEP warnings when CASTEP not installed
from tkinter import *
from tkinter import filedialog # used for importing files
from termcolor import cprint, colored
from pyfiglet import figlet_format
from colorama import init, Fore, Style # used for coloured output in terminal
init(convert=True, strip=not sys.stdout.isatty()) # initialises colour conversion in terminal


#########################################################################################
### DISPLAY SCRIPT NAME
cprint(figlet_format(' Cell Builder ', font='slant'), 'cyan')


#########################################################################################
### FILE BROWSER FOR USER TO SELECT A .CIF (or other ASE-supported file)
cprint('PICK A .CIF:', 'red')
input('A pop-up window will open where you can search for your file. '
	'Click Enter to continue.')
print()

# use tkinter to create a file browser window
file_not_selected = True
while file_not_selected:
    root = Tk()
    root.title('File Picker')
    root.geometry('500x10') # root window, small because not necessary
    filename = filedialog.askopenfilename(title="File Picker")
    root.destroy()
    # check if a file was selected; file selection is not forced
    if len(filename) != 0:
        file_not_selected = False
    else:
        cprint('File not selected. Try again.', 'red')
        print()

# print the path to the file
print(colored('File chosen: ', 'cyan'), colored(filename, 'yellow'))
print()

#########################################################################################
### ASK FOR DETAILS ABOUT THE SYSTEM
cprint('GIVE COIL ORIENTATION:', 'red')
print('Set the crystallographic axis aligned with the direction of the coil (a/b/c).')
print()

abc_dict = {'a':0, 'b':1, 'c':2} # interprets coil axis
pq_dict = {0:[1,2], 1:[0,2], 2:[0,1]}

# ask for coil axis
coil_axis_not_selected = True
while coil_axis_not_selected:
	cprint('Coil axis = ', 'cyan', end='')
	coil_axis_input = input(Fore.YELLOW)
	# interpret user input
	if coil_axis_input in abc_dict.keys():
		# define coil, p and q axes
		coil_axis = abc_dict[coil_axis_input]
		p_axis = pq_dict[coil_axis][0]
		q_axis = pq_dict[coil_axis][1]

		coil_axis_not_selected = False
	else:
		cprint('Incorrect coil axis! Try again.', 'red')
		print()

print(Style.RESET_ALL)

#########################################################################################
# READ THE FILE AND CREATE VARIABLES
with warnings.catch_warnings():
	warnings.simplefilter('ignore')
	atoms = ase.io.read(filename) # read structural information into ASE Atoms object

# get useful information from atoms object
symbols = atoms.get_chemical_symbols() # chemical symbols as a list
cell_lengths = atoms.get_cell().lengths()
cell_angles = atoms.get_cell().angles()

cell_angles[coil_axis] = 90

new_cell_params = np.concatenate((cell_lengths, cell_angles))

atoms.set_cell(new_cell_params)

atoms.center()


cell_matrix = atoms.get_cell()
cart_positions = atoms.get_positions() # cartesian position vectors as a list
frac_positions = atoms.get_scaled_positions() # fractional positions as a list

# angle phi between the normal to the pc plane and the q vector;
# this is used to calculate the new unit cell q axis length from the
# required distance between coils
phi = np.arccos((
		(np.dot(np.cross(cell_matrix[coil_axis], cell_matrix[p_axis]),
				cell_matrix[q_axis]))
        /
		(np.linalg.norm(np.cross(cell_matrix[coil_axis], cell_matrix[p_axis]))
		* np.linalg.norm(cell_matrix[q_axis]))))

# angle xi between the normal to the qc plane and the p vector;
# this is used to calculate the new unit cell p axis length from the
# required distance between coils
xi = np.arccos((
		(np.dot(np.cross(cell_matrix[coil_axis], cell_matrix[q_axis]),
				cell_matrix[p_axis]))
        /
		(np.linalg.norm(np.cross(cell_matrix[coil_axis], cell_matrix[q_axis]))
		* np.linalg.norm(cell_matrix[p_axis]))))

#########################################################################################
# GIVE INFORMATION ABOUT COIL WIDTH

frac_positions_to_fit_list = list()
cart_positions_to_check = list()
for i in range(0, len(symbols)):
    if symbols[i] in {'C', 'N', 'O'}:
        frac_positions_to_fit_list.append(frac_positions[i])
        cart_positions_to_check.append(cart_positions[i])

frac_positions_to_fit = np.array(frac_positions_to_fit_list)

avg_coil_frac = np.mean(frac_positions_to_fit, axis=0)
avg_coil_cart = cell_matrix.cartesian_positions(avg_coil_frac)

coil_line = Line(point=Point(avg_coil_cart), direction=cell_matrix[coil_axis])

coil_radius = 0
for i in range(0, len(cart_positions_to_check)):
	distance = coil_line.distance_point(cart_positions_to_check[i])
	if distance > coil_radius:
		coil_radius = distance

cprint('At its widest, the coil has a radius of: ', 'cyan', end='')
cprint("{:.6f}".format(coil_radius), 'yellow', end='')
cprint(' angstroms.', 'cyan')
print()

#########################################################################################
# ASK FOR DISTANCES TO SET
cprint('SET DISTANCE:', 'red')
print('Input the desired distances between coils.')
print()

cprint('Coil-coil distances = ', 'cyan', end='')
cc_distances_input = input(Fore.YELLOW)

cc_distances = tuple(float(item) for item in cc_distances_input.split(','))

print()

#########################################################################################
# ASK WHETHER TO ADD COUNTER-IONS AT WALLS
cprint('ADD COUNTER-IONS?', 'red')
print('Do you want to add counter-ions at the walls of the cell? (y/n)')
print()

add_ions = None
while add_ions == None:
	cprint('Counter-ions = ', 'cyan', end='')
	ions_input = input(Fore.YELLOW)
	# interpret user input
	if ions_input in {'y', 'Y', 'yes', 'Yes', 'YES'}:
		add_ions = True
	elif ions_input in {'n', 'N', 'no', 'No', 'NO'}:
		add_ions = False
	else:
		cprint('Incorrect input! Try again.', 'red')
		print()

print(Style.RESET_ALL)



#########################################################################################
### WRITE FILES FOR EACH PROVIDED DISTANCE

counter = 1

for cc_distance in cc_distances:

	p_axis_length = cc_distance / abs(np.cos(xi))
	q_axis_length = cc_distance / abs(np.cos(phi))

	new_cell_params[p_axis] = p_axis_length
	new_cell_params[q_axis] = q_axis_length

	new_atoms = Atoms(symbols=symbols, cell=new_cell_params, positions=cart_positions)
	new_atoms.center()
	new_cell_matrix = new_atoms.get_cell()
	

	if add_ions:
		counter_ion_1 = np.zeros(3) # in the pc plane
		counter_ion_1[q_axis] = 0.00
		counter_ion_1[p_axis] = 0.50
		counter_ion_1[coil_axis] = 0.50

		counter_ion_2 = np.zeros(3) # in the qc plane
		counter_ion_2[q_axis] = 0.50
		counter_ion_2[p_axis] = 0.00
		counter_ion_2[coil_axis] = 0.50

		frac_positions_to_write = [counter_ion_1, counter_ion_2]
		symbols_to_write = ['Cl', 'Cl']
	else:
		frac_positions_to_write = list()
		symbols_to_write = list()


	for position in new_atoms.get_scaled_positions():
			frac_positions_to_write.append(position)

	for symbol in new_atoms.get_chemical_symbols():
			symbols_to_write.append(symbol)

	if add_ions:
		new_filename = f'{filename.split("_input.cell")[0]}_ci_{cc_distance/2}.cell'
	else:
		new_filename = f'{filename.split("_input.cell")[0]}_charged_{cc_distance}.cell'
	# write new cell file    
	with open(new_filename, 'w') as cell_file:
		# write cell vectors
		cell_file.write('%BLOCK lattice_cart')
		cell_file.write('\n')
		cell_file.write("{:.6f}".format(new_cell_matrix[0,0]) + '\t'
						+ "{:.6f}".format(new_cell_matrix[0,1]) + '\t'
						+ "{:.6f}".format(new_cell_matrix[0,2]))
		cell_file.write('\n')
		cell_file.write("{:.6f}".format(new_cell_matrix[1,0]) + '\t'
						+ "{:.6f}".format(new_cell_matrix[1,1]) + '\t'
						+ "{:.6f}".format(new_cell_matrix[1,2]))
		cell_file.write('\n')
		cell_file.write("{:.6f}".format(new_cell_matrix[2,0]) +'\t'
						+ "{:.6f}".format(new_cell_matrix[2,1]) + '\t'
						+ "{:.6f}".format(new_cell_matrix[2,2]))
		cell_file.write('\n')
		cell_file.write('%ENDBLOCK lattice_cart')
		cell_file.write('\n')
		cell_file.write('\n')
		# write fix all cell boolean
		cell_file.write('FIX_ALL_CELL : TRUE')
		cell_file.write('\n')
		cell_file.write('\n')
		cell_file.write('\n')
		# write atomic coordinates
		cell_file.write('%BLOCK positions_frac')
		cell_file.write('\n')
		for i in range(0, len(frac_positions_to_write)):
			cell_file.write(symbols_to_write[i]  + '\t' +
						"{:.10f}".format(frac_positions_to_write[i][0]) + '\t' +
						"{:.10f}".format(frac_positions_to_write[i][1]) + '\t' +
						"{:.10f}".format(frac_positions_to_write[i][2]))
			cell_file.write('\n')
		cell_file.write('%ENDBLOCK positions_frac')
		cell_file.write('\n')
		cell_file.write('\n')
       	# write species pot
		cell_file.write('%BLOCK species_pot')
		cell_file.write('\n')
		cell_file.write('%ENDBLOCK species_pot')
		cell_file.write('\n')
		cell_file.write('\n')
       	# write symmetry generate command
		cell_file.write('SYMMETRY_GENERATE')
		cell_file.write('\n')
		cell_file.write('\n')
       	# write k-point mp spacing
		cell_file.write('KPOINT_MP_SPACING : 0.05')

counter += 1
  


print()
print()
cprint('Done. Exiting.', 'green')

print(Style.RESET_ALL)


