
import streamlit as st
from stat_thermo import *

def add_molecule(sign):
    name = st.text_input('Name/Formula', help='Name of the molecule')
    s_c = float(st.text_input('stoichiometric number', help='Reaction stoichiometric number (positive value)'))
    P = float(st.text_input('P: Pressure [bar]'))
    m = float(st.text_input('m: molecular mass [dalton]'))
    B = float(st.text_input('B: rotational constant [cm-1]'))
    A = float(st.text_input('A: other rotational constant [cm-1]'), 0)
    C = float(st.text_input('C: other rotational constant [cm-1]'), 0)
    linearity_dict = {'Linear': True, 'Non-linear': False, 'an Atom': 'Atom'}
    linear = linearity_dict[ st.selectbox('The molecule is', ('Linear', 'Non-linear', 'an Atom')) ]
    o = float(st.text_input('o: symmetry number', help='The symmetry number a of a molecule is the order of the finite rotational sub-group of the point group of the molecule'))
    n_mod_raw = st.text_input('List of the vibrational mode frequencies [cm-1] separated by a comma').split(',')
    n_mod = [float(i) for i in n_mod_raw]
    n_deg_raw = st.text_input('List of the vibrational mode degeneracies separated by a comma').split(',')
    n_deg = [float(i) for i in n_deg_raw]
    En_raw = st.text_input('List of the electron level Energies [cm-1] separated by a comma').split(',')
    En = [float(i) for i in En_raw]
    gn_list_elec_raw = st.text_input('List of the electron levels degeneracies separated by a comma').split(',')
    gn_list_elec = [float(i) for i in gn_list_elec_raw]
    spin_list_raw = st.text_input('List of nuclear spins separated by a comma').split(',')
    spin_list = [float(i) for i in spin_list_raw]
    load = st.button('submit molecule')
    while not load:
        if load:
            data[name]["param"] = [P, m, B, o, linear, n_mod, n_deg, gn_list_elec, En, spin_list, A, C]




if __name__ == "__main__":
    st.write('''
    # Statistical Thermodynamics Reactor!
    ''')

    st.write('''
    Calculates the thermodynamics state functions for a reaction!
    For more options and information, check out the [GitHub repository](https://github.com/lmonari5/stat_thermo_reactor.git)
    ''')

    dU0 = float(st.text_input('ΔU0 [Joule]', help='ΔU at 0 Kelvin'))
    T = float(st.text_input('Temperature [Kelvin]'))

    st.write('''
    ## Reactants:
    ''')
    st.button('Add', on_click=add_molecule())

    print(data)


