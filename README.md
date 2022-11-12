# Statistical Thermodynamics Reactor
A python script to calculate reaction data with Statistical Thermodynamics (ΔG, ΔS, ΔH, ΔU, reaction constant); from the lecture of Roberto Marquardt.

## Try the Web App:

Try the web app based on Streamlit: it's free, no need to install anything!

[Stat-Therm-Reactor Web App](https://lmonari5-stat-thermo-reactor-streamlit-app-2c05u1.streamlitapp.com)

### Required package:
- pandas
- matplotlib

### Remarks:
- The script works for one reaction only. 
- The input parameters are stored in a dictionary called 'data'. The reaction state functions are stored in a dictionary called 'results'.
- The script is meant to be used from python interpreter. Make sure the script is located in the working directory or
added to the pythonpath
- To convert cm-1 to Joule: multiply for hc

  
### Required parameters summary:
- U0_value: Internal Energy of the reaction at 0 kelvin at standard pressure [Joule/mol]
- s_c: stoichiometric number, positive for products and negative for reactants
- P: pressure [bar]
- m: molecular mass [dalton]
- B: rotational constant [cm-1]
- o: symmetry number of the compound (The symmetry number of a molecule is the order of the finite rotational sub-group of the point group of the molecule.)
- linear: the molecule is linear ('True'), not linear ('False'), or an atom ('Atom')
- n_mod: list of the vibrational mode frequencies [cm-1]
- n_deg: list of the vibrational mode degeneracies
- gn_list_elec: list of electron levels degeneracies
- En: list of electron level Energies [cm-1]
- spin_list: list of the nuclear spins
- A: other eventual rotational constant [cm-1]
- C: other eventual rotational constant [cm-1]
- T: temperature [Kelvin]


## Usage in the python interpreter:

Import the functions: \
`from stat_thermo import *`

Initialize a dictionary with the value of the Internal Energy at 0 Kelvin (`U0_value`), by typing

```
data["U0"]= U0_value
``` 

Initialize a new compound (`compound_name`) and set the stoichiometric number:

```
data["compound_name"]=dict()
data["compound_name"]["s_c"]= stechometric number
```

Start setting the parameters for each reactant and product. Do not change the keyword 'param'. ***Use this order:***

`data['compound_name']["param"]= [P, m, B, o, linear (True, False or 'Atom'), n_mod, n_deg, gn_list_elec, En, spin_list, A, C]`

## Functions

#### Low-level functions:
- Partition functions:
    - qt(T, P, m)
    - qro(T, B, linear, o, A=0, C=0)
    - qv(T, n_mod, n_deg)
    - qe(T, gn_list_elec, En)
    - qn(spin_list)
    - Q(name, T, P, m, B, o, linear, n_mod, n_deg, gn_list_elec, En, spin_list, A=0, C=0)
- Internal Energy:
  - Ut(T)
  - Ur(T, linear)
  - Uv(T, n_mod, n_deg)
  - Ue(T, qe, gn_list_elec, En)
  - U(compound_name, T, n_mod, n_deg, linear, qe, gn_list_elec, En)
- Entropy:
  - St(qt)
  - S(T, q, U)
  - S_calc(compound_name, T)

#### Easy-to-use functions:
- Q_fast(compound_name, T)
- U_fast(name, T)
- S_fast(name, T)
- k_fast(T)

#### Plotting the state functions:
- plot(Ti, Tf, step, save=False)

#### Print the results on screen:
- visual()

#### Save the results:
- save(filename='stat_thermo_results')

## Example:
Decomposition of formaldehyde  into carbon monoxide and dihydrogen:

1 H2CO ⇌ 1 H2 + 1 CO
```
# import the functions

from stat_thermo import *

# set up the initial parameters

data["U0"]= -1700
data["H2"]=dict()
data["H2"]["s_c"]= 1
data['H2']['param']=[1, 2.01588, 60.89, 2, True,[4280], [1], [1], [0], [1/2, 1/2], 0, 0]
data["CO"]=dict()
data["CO"]["s_c"]= 1
data["CO"]["param"]=[1, 28.0140, 1.930, 1, True, [2156], [1], [1], [0], [0, 0], 0, 0]
data["H2CO"]=dict()
data["H2CO"]["s_c"]= -1
data["H2CO"]["param"]=[1, 30.02628, 9.405, 2, False, [2843, 2766, 1746, 1501, 1247, 1164],[1,1,1,1,1,1], [1], [0], [1/2, 1/2, 0, 0], 1.295, 1.134 ]

# Calculate the results

Q_fast("H2",300) #calculate all q for H2
data["H2"][300]["qt"] #to see the values of qt
U_fast("H2",300)
S_fast("H2",300)
k_fast(300)
k_fast(100)
visual()
plot(100,300,20, save= "example")
visual()
save("reaction_data")
```

## Help

In the python interpreter type:
```
import stat_thermo
help(stat_thermo)
```
