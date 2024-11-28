
# README
## zadanie1.py
### Usage
`zadanie1.py [-h] [--name STRUCTURE_NAME]`

#### Running with Default Values
To run the script using default values (structure_name = 4YWO), use:
```bash
python zadanie1.py 
```
---
## zadanie2.py


### Usage
`zadanie2.py [-h] [--name STRUCTURE_NAME] [--function FUNCTION USED TO CALCULATE ANGLES]`

#### Running with Default Values
To run the script using default values (structure_name = 4YWO, function = calc_dihedral), use:
```bash
python zadanie2.py 
```
#### Functions available for calculating angles
```
.calc_dihedral - [--function calc_dihedral]
.get_phi_psi_list - [--function get_phi_psi_list]
```
---

### Examples Commend

```bash
python .\zadanie2.py --name 2hhb --function get_phi_psi_list
```

## main.py

### Usage
```bash
python main.py [-h] [--load_from_disc LOAD_FROM_DISC] [--output OUTPUT]
               [--oid OID]
               structure_id
```
### Example use
```bash
python main.py 430d.pdb --output 430d_cg.pdb
```
file used with 430d.pdb is /examples/430d.pdb
output file is /examples/430d_cg.pdb

## coarse_grain_rebuiler
### Usage
```bash
python coarse_grain_rebuilder.py [-h] [--output_file OUTPUT_FILE] [--oid OID]
                                 [--templates TEMPLATES]
                                 input
```
### example use
python .\coarse_grain_rebuilder.py .\examples\430d_cg.pdb --output_file .\examples\430d_r.pdb --oid 430dr

input and output files are in folder examples