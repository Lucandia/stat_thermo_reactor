# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 13:27:58 2021

@author: lmonari5
"""
import math
import matplotlib.pyplot as plt
import pandas as pd

# useful constant:
kb = 1.38064852e-23
h = 6.62607004e-34
c = 2.99792458e10  # IN CENTIMETER
Na = 6.0221409e23
c2 = h * c / kb
# convert cm-1 to joule: *hc
data = dict()
global results
results = dict()

'''
Initialize a dictionary of reactant and product, by typing
data = dict()
data["U0"]= U0_value

data["compound_name"]=dict()

Set the stechiometric number by tiping:
data["compound_name"]["s_c"]= stechometric number

Insert the paramter in this order (important), by tiping:
data['compound_name']["param"]= [P, m, B, o, linear (True, False or 'Atom'), n_mod, n_deg, gn_list_elec, En, spin_list, A, C]

example:
data = dict()
data["U0"]= -1700

data["H2"]=dict()
data["H2"]["s_c"]= 1
data['H2']["param"]=[1, 2.01588, 60.89, 2, True,[4280], [1], [1], [0], [1/2, 1/2], 0, 0]

Q_fast('H2',300)

S_fast('H2', 300) #to see just the data for H2

data["CO"]=dict()
data["CO"]["s_c"]= 1
data['CO']["param"]=[1, 28.0140, 1.930, 1, True, [2156], [1], [1], [0], [0, 0], 0, 0]

data["H2CO"]=dict()
data["H2CO"]["s_c"]= -1
data['H2CO']["param"]=[1, 30.02628, 9.405, 2, False, [2843, 2766, 1746, 1501, 1247, 1164],[1,1,1,1,1,1], [1], [0], [1/2, 1/2, 0, 0], 1.295, 1.134 ]

k_fast(300) #calculate everything on his own
k_fast(100)
visual()
plot(100,300,1, save='last_exercise')
save('reaction_data')
'''

print("IMPORTANT: USE THIS PROGRAM JUST FOR 1 REACTION\n"
      'Initialize a dictionary of reactant and product, by typing\n' \
      'data = dict()\n' \
      'data["U0"]= U0_value\n' \
      'data["compound_name"]=dict()\n\n' \
      'Set the stechiometric number by tiping:\n' \
      'data["compound_name"]["s_c"]= stechometri number\n\n' \
      'Insert the paramter in this order, by tiping:\n' \
      "data['compound_name']['param']= [Pressure (bar), mass (dalton), B (rot_const cm-1)," \
      'symmetry_number, linear (True, False or "Atom"), list_vibrat_freq (cm-1), list_vibrat_degeneracy,' \
      'list_elec_degeneracy, list_elec_energy (cm-1), spin_list_nuclei, A (rot_const), C (rot_const)]\n\n' \
      'With "Atom" option qrot and qvib are set equal to 1 (they do not have a meaning for atoms)'
      'To know the value of a calculated variable type:\n' \
      'data["compound_name"][Temperature]["variable"]\n\n' \
      'Function of interest:\n' \
      "Q_fast('compound_name',T) #calculate all Q for a compound from scratch \n" \
      "S_fast('compound_name',T) #calculate all S for a compount from scratch \n" \
      "k_fast(T) #calculate everything directly\n" \
      "plot(initial_T,final_T, steps, save='filename') #calculate each K and plot the results for every variable\n" \
      "visual() visualize all the results in a table\n" \
      "save(filename) #save the data in txt and excel files\n" \
      "for other functions or details look into the code\n"

      '\nexample:\n' \
      'data["U0"]= -1700\n' \
      'data["H2"]=dict()\n' \
      'data["H2"]["s_c"]= 1\n' \
      "data['H2']['param']=[1, 2.01588, 60.89, 2, True,[4280], [1], [1], [0], [1/2, 1/2], 0, 0]\n" \
      'data["CO"]=dict()\n' \
      'data["CO"]["s_c"]= 1\n' \
      'data["CO"]["param"]=[1, 28.0140, 1.930, 1, True, [2156], [1], [1], [0], [0, 0], 0, 0]\n' \
      'data["H2CO"]=dict()\n' \
      'data["H2CO"]["s_c"]= -1\n' \
      'data["H2CO"]["param"]=[1, 30.02628, 9.405, 2, False, [2843, 2766, 1746, 1501, 1247, 1164],[1,1,1,1,1,1], [1], [0], [1/2, 1/2, 0, 0], 1.295, 1.134 ]\n' \
      '\n#Then you can calculate the value of interest\n' \
      'Q_fast("H2",300)\n #calculate all q for H2\n' \
      'data["H2"][300]["qt"]=dict()\n #to see the values of qt\n' \
      'U_fast("H2",300)\n' \
      'S_fast("H2",300)\n' \
      'k_fast(300)\n' \
      'k_fast(100)\n' \
      'visual()\n' \
      'plot(100,300,1, save= "last_excercise")\n' \
      'visual()\n' \
      'save("reaction_data")\n')


# translational partition function as a function of T,P,N, using unit measure of
# SI and m in dalton
# ***qt divided by N*** qt_norm
def qt(T, P, m):
    """> This function calculates the canonical translational partition function of a compound at a given temperature and pressure

    Parameters
    ----------
    T
    	Temperature in Kelvin
    P
    	Pressure in bar
    m
    	Mass in Dalton

    Returns
    -------
    	The translational partition function

    """
    return 0.02594678522 * m ** (3 / 2) * T ** (5 / 2) / P


def qr(T, B, gn_even=1, gn_odd=1, threshold=1e-6):
    """> This function calculates the canonical rotational partition function (without approximation) of a compound at a given
    temperature, rotational constant and levels degeneracy.

    Parameters
    ----------
    T
    	Temperature in Kelvin
    B
    	Rotational constant in cm-1
    gn_even, optional
    	Degeneracy of even-parity levels
    gn_odd, optional
    	Degeneracy of odd-parity levels
    threshold, optional (1e-6 by default)
    	The function stops when the contribution of new levels is below the threshold

    Returns
    -------
    	The rotational partition function

    """
    qrtot = float()
    qr = 1
    j = 0
    while qr > threshold:  # precision on the threshold
        if j % 2 == 0:
            qr = gn_even * (2 * j + 1) * math.exp(-h * c * B * j * (j + 1) / (kb * T))
        else:
            qr = gn_odd * (2 * j + 1) * math.exp(-h * c * B * j * (j + 1) / (kb * T))
        qrtot += qr
        j += 1
    return qrtot


def qro(T, B, linear, o, A=0, C=0):
    """> This function calculates the canonical rotational partition function from the symmetry number (without approximation)
    of a compound at a given temperature, rotational constants and symmetry number

    Parameters
    ----------
    T
    	Temperature in Kelvin
    B
    	Main rotational constant in cm-1
    linear
    	Boolean value that determines whether the system is linear or not
    o
    	Symmetry number of the molecule (The symmetry number a of a molecule is the order of the finite rotational sub-group of the point group of the molecule)
    A, optional
    	Rotational constant in cm-1
    C, optional
    	Rotational constant in cm-1

    Returns
    -------
    	The rotational partition function

    """
    if linear:
        qr = kb * T / (h * c * B * o)
    else:
        qr = math.sqrt(math.pi / (A * B * C)) * (kb * T / (h * c)) ** (3 / 2) / o
    return qr


def qv(T, n_mod, n_deg):
    """> This function calculates the canonical vibrational partition function of a compound at a given temperature,
    considering the vibrational modes frequencies and degeneracies

    Parameters
    ----------
    T
    	Temperature in Kelvin
    n_mod
    	List of vibrational modes frequencies
    n_deg
    	List of vibrational modes degeneracies

    Returns
    -------
    	The vibrational partition function

    """
    if type(n_mod) != type(list()):
        return print('error, submit a list of frequency for modes')
    qv = 1
    for index, vs in enumerate(n_mod):
        qs = (1 / (1 - math.exp(-vs * c2 / (T)))) ** n_deg[index]
        qv *= qs
    return qv


def qe(T, gn_list, En):
    """> This function calculates the canonical electronic partition function of a compound at a given temperature,
    considering the electronic modes frequencies and degeneracies

    Parameters
    ----------
    T
    	Temperature in Kelvin
    gn_list
    	List of electronic levels degeneracies
    En
    	List of electronic levels energies

    Returns
    -------
    	The electronic partition function

    """
    qe = 0
    for index, gn in enumerate(gn_list):
        qe += gn * math.exp(-En[index] * c2 / T)
    return qe


def qn(si_list):
    """> This function calculates the canonical nuclear partition function of a compound at a given temperature,
    considering the spin quantum number for each nucleus

    Parameters
    ----------
    si_list
    	List of the spin quantum number for each nucleus

    Returns
    -------
    	The nuclear partition function

    """
    qn = 1
    for si in si_list:
        qn *= (2 * si + 1)
    return qn


def Ezpv(vs_list, ds_list=1):
    """> This function takes a list of vibration frequencies and a list of degeneracies and returns the zero point energy

    Parameters
    ----------
    vs_list
    	list of vibration frequencies
    ds_list, optional
    	degeneracy of each vibrational mode

    Returns
    -------
    	The zero point vibrational energy

    """
    E = 0
    if ds_list == 1:
        print('all vibration degeneracy are equal to 1')
        ds_list = len(vs_list) * '1'
    for index, vs in enumerate(vs_list):
        E += 1 / 2 * h * c * vs * float(ds_list[index])
    return E


# all partition function in once:
# print('name,T, P, m, B, o, linear, n_mod, n_deg, gn_list_elec, En, spin_list, A, C')
def Q(name, T, P, m, B, o, linear, n_mod, n_deg, gn_list_elec, En, spin_list, A=0, C=0):
    """> This function calculates all the partition function for a given molecule

    Parameters
    ----------
    name
    	Name of the molecule
    T
    	Temperature in Kelvin
    P
    	Pressure in bar
    m
    	Mass of the molecule in dalton
    B
    	Rotational constant
    o
    	Symmetry number of the molecule
    linear
    	True if the molecule is linear, False if it is not, Atom if it is an atom
    n_mod
    	List of vibrational modes frequencies
    n_deg
    	List of vibrational modes degeneracies
    gn_list_elec
    	List of electronic levels degeneracies
    En
    	List of electronic levels energies
    spin_list
    	List of the spin quantum number for each nucleus
    A, optional
    	Rotational constant
    C, optional
    	Rotational constant

    Returns
    -------
    	The total partition function

    """
    qti = qt(T, P, m)
    if linear == 'Atom':
        qri = 1
        qvi = 1
    else:
        qri = qro(T, B, linear, o, A, C)
        qvi = qv(T, n_mod, n_deg)
    qei = qe(T, gn_list_elec, En)
    qni = qn(spin_list)
    Q = qti * qri * qvi * qei * qni
    data[name][T] = dict()
    data[name][T]['qt'] = qti
    data[name][T]['qr'] = qri
    data[name][T]['qv'] = qvi
    data[name][T]['qe'] = qei
    data[name][T]['qn'] = qni
    data[name][T]['Q'] = Q
    return Q


def Ut(T):
    """> This function calculates the translational internal energy of the system.

    Parameters
    ----------
    T
    	Temperature in Kelvin

    Returns
    -------
    	The translational internal energy

    """
    return 3 / 2 * kb * T * Na


def Ur(T, linear):
    """> This function calculates the rotational internal energy of the system.

    Parameters
    ----------
    T
    	Temperature in Kelvin
    linear
    	Boolean value that determines whether the system is linear or not

    Returns
    -------
    	The rotational internal energy

    """
    if linear:
        U = kb * T
    else:
        U = 3 / 2 * kb * T
    return U * Na


def Uv(T, n_mod, n_deg):
    """> This function calculates the vibrational internal energy of the system.

    Parameters
    ----------
    T
    	Temperature in Kelvin
    n_mod
    	List of vibrational modes frequencies
    n_deg
    	List of vibrational modes degeneracies

    Returns
    -------
    	The vibrational internal energy

    """
    if type(n_mod) != type(list()):
        return print('error, submit a list of frequency for modes')
    Uv = 0
    for index, vs in enumerate(n_mod):
        Us = n_deg[index] * h * vs * c / (math.exp(c2 * vs / T) - 1)
        Uv += Us
    return Uv * Na


def Ue(T, qe, gn_list_elec, En):
    """> This function calculates the electronic internal energy of the system.

    Parameters
    ----------
    T
    	Temperature in Kelvin
    qe
        Electronic partition function
    gn_list_elec
    	List of electronic levels degeneracies
    En
    	List of electronic levels energies

    Returns
    -------
    	The electronic internal energy

    """
    Ue = 0
    for index, gn in enumerate(gn_list_elec):
        Ue += gn * En[index] * h * c * math.exp(-En[index] * c2 / T) / qe
    return Ue * Na


def U(name, T, n_mod, n_deg, linear, qe, gn_list_elec, En):
    """> This function calculates the translational internal energy of the system.

    Parameters
    ----------
    name
    	Name of the molecule
    T
    	Temperature in Kelvin
    n_mod
    	List of vibrational modes frequencies
    n_deg
    	List of vibrational modes degeneracies
    linear
    	Boolean value that determines whether the system is linear or not
    qe
        Electronic partition function
    gn_list_elec
    	List of electronic levels degeneracies
    En
    	List of electronic levels energies

    Returns
    -------
    	The total internal energy

    """
    Uti = Ut(T)
    if linear == 'Atom':
        Uri = 0
        Uvi = 0
    else:
        Uri = Ur(T, linear)
        Uvi = Uv(T, n_mod, n_deg, )
    Uei = Ue(T, qe, gn_list_elec, En)
    Utot = Uti + Uri + Uvi + Uei
    data[name][T]['Ut'] = Uti
    data[name][T]['Ur'] = Uri
    data[name][T]['Uv'] = Uvi
    data[name][T]['Ue'] = Uei
    data[name][T]['Utot'] = Utot
    return Utot


def S(T, q, U):
    """> This function calculates the entropy of the system given a partition function.

    Parameters
    ----------
    T
    	Temperature in Kelvin
    q
    	A partition function
    U
    	Internal energy

    Returns
    -------
    	The entropy of the system according to the given partition function

    """
    return Na * kb * math.log(q) + (U) / T


def St(qt):
    """> This function calculates the translational entropy of the system according to the Sackur-Tetrode equation.

    Parameters
    ----------
    qt
    	The translational partition function

    Returns
    -------
    	The translational entropy

    """
    return Na * kb * math.log(qt) + 5 / 2 * Na * kb


def S_calc(name, T):
    """> This function calculates the entropy of a molecule stored in the data dictionary and store them in the
        dictionary. Before running this function, you have to calculate the internal energies and partition functions

    Parameters
    ----------
    name
    	The name of the compound stored in data
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    Sti = St(data[name][T]['qt'])
    Sr = S(T, data[name][T]['qr'], data[name][T]['Ur'])
    Sv = S(T, data[name][T]['qv'], data[name][T]['Uv'])
    Se = S(T, data[name][T]['qe'], data[name][T]['Ue'])
    Stot = Sti + Sr + Sv + Se
    data[name][T]['St'] = Sti
    data[name][T]['Sr'] = Sr
    data[name][T]['Sv'] = Sv
    data[name][T]['Se'] = Se
    data[name][T]['Stot'] = Stot
    return


# list of Q of compounds and list of stoichiometric number
def k(T):
    """> This function calculates the reaction constant of all the species in the 'data' dictionary and store it in the
        'results' dictionary. The constant is calculated from the partition functions

    Parameters
    ----------
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    results[T] = dict()
    k = math.exp(data['U0'] / (kb * Na * T))
    # print('Assuming Na molecules')
    for compound in data:
        if compound == 'U0': continue
        stec_coeff = data[compound]['s_c']
        Q = data[compound][T]['Q']
        k *= Q ** stec_coeff
        results[T]['k'] = k
    return


def Q_fast(name, T):
    """> This function calculates the partition functions of a molecule stored in the data dictionary and store them in the
            dictionary

    Parameters
    ----------
    name
    	The name of the compound stored in data
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    print('\n' + name)
    P = data[name]["param"][0]
    m = data[name]["param"][1]
    B = data[name]["param"][2]
    o = data[name]["param"][3]
    linear = data[name]["param"][4]
    n_mod = data[name]["param"][5]
    n_deg = data[name]["param"][6]
    gne_list = data[name]["param"][7]
    En = data[name]["param"][8]
    spin_list = data[name]["param"][9]
    A = data[name]["param"][10]
    C = data[name]["param"][11]
    qti = qt(T, P, m)
    if linear == 'Atom':
        qri = 1
        qvi = 1
    else:
        qri = qro(T, B, linear, o, A, C)
        qvi = qv(T, n_mod, n_deg)
    qei = qe(T, gne_list, En)
    qni = qn(spin_list)
    Q = qti * qri * qvi * qei * qni
    print('qt', qti)
    print('qrot', qri)
    print('qv', qvi)
    print('qe', qei)
    print('qn', qni)
    print('Q:', Q)
    data[name][T] = dict()
    data[name][T]['qt'] = qti
    data[name][T]['qr'] = qri
    data[name][T]['qv'] = qvi
    data[name][T]['qe'] = qei
    data[name][T]['qn'] = qni
    data[name][T]['Q'] = Q
    return Q


def U_fast(name, T):
    """> This function calculates the internal energies of a molecule stored in the data dictionary and store them in the
            dictionary. This function calculates on his own the partition functions

    Parameters
    ----------
    name
    	The name of the compound stored in data
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    Q_fast(name, T)
    U(name, T, data[name]["param"][5], data[name]["param"][6], data[name]["param"][4],
      data[name][T]['qe'], data[name]["param"][7], data[name]["param"][8])
    print('Ut', data[name][T]['Ut'])
    print('Ur', data[name][T]['Ur'])
    print('Uv', data[name][T]['Uv'])
    print('Ue', data[name][T]['Ue'])
    print('Utot', data[name][T]['Utot'])
    return


def S_fast(name, T):
    """> This function calculates the entropy of a molecule stored in the data dictionary and store them in the
            dictionary. This function calculates on his own the partition functions and internal energies

    Parameters
    ----------
    name
    	The name of the compound stored in data
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    U_fast(name, T)
    S_calc(name, T)
    print('St', data[name][T]['St'])
    print('Sr', data[name][T]['Sr'])
    print('Sv', data[name][T]['Sv'])
    print('Se', data[name][T]['Se'])
    print('Stot', data[name][T]['Stot'])
    return


def k_fast(T):
    """> This function calculates the reaction constant of all the species in the 'data' dictionary and store it in the
        'results' dictionary. This functions autonomously calculates all the Thermodynamics state functions (U,S,H,G),
        store them in the 'results' dictionary and show that the dG obtained by the state functions or by the reaction
        constants should match.

    Parameters
    ----------
    T
        The temperature in Kelvin

    Returns
    -------
    	None

    """
    results[T] = dict()
    print('\n' + str(T) + ' K:')
    k = math.exp(-data['U0'] / (kb * Na * T))
    deltaU = float()
    delta_s_c = float()
    deltaS = float()
    for compound in data:
        if compound == 'U0': continue
        S_fast(compound, T)
        k *= data[compound][T]['Q'] ** data[compound]['s_c']
        deltaU += data[compound][T]["Utot"] * data[compound]['s_c']
        delta_s_c += data[compound]['s_c']
        deltaS += data[compound][T]["Stot"] * data[compound]['s_c']
        for var in data[compound][T]:
            if 'U' in var and 'tot' not in var or 'S' in var and 'tot' not in var:
                print(var, data[compound][T][var])
                if int(data[compound][T][var]) == 0:
                    results[T][var] = results[T].get(var, 0) + data[compound][T][var]
                    continue
                results[T][var] = results[T].get(var, 0) + data[compound][T][var] * data[compound]['s_c']
    deltaU += +data['U0']
    deltaH = deltaU + Na * kb * T * delta_s_c
    deltaG_thermo = deltaH - T * deltaS
    deltaG_k = -Na * kb * T * math.log(k)
    print('\nreaction data\n')
    print('k', k)
    print("\u0394U", deltaU)
    print("\u0394S", deltaS)
    print("\u0394H", deltaH)
    print("\u0394G thermo", deltaG_thermo)
    print("\u0394G from k", deltaG_k)
    print('the two \u0394G should be equal')
    results[T]
    results[T]['k'] = k
    results[T]['\u0394U'] = deltaU
    results[T]['\u0394S'] = deltaS
    results[T]['\u0394H'] = deltaH
    results[T]['\u0394G_thermo'] = deltaG_thermo
    results[T]['\u0394G_k'] = deltaG_k
    return


def plot(Ti, Tf, step, save=False):
    """>  This function calculates and plot the state functions from Temperature Ti to Temperature Tf with a dT equal to
        the step.

    Parameters
    ----------
    Ti
        The initial temperature in Kelvin
    Tf
        The final temperature in Kelvin
    step
        The dT in Kelvin
    save, optional
        Boolean value, wether to save the plots or not.

    Returns
    -------
    	None

    """
    T = Ti
    T_list = list()
    while T <= Tf:
        k_fast(T)
        T_list.append(T)
        T += step
    for delta in ['\u0394U', '\u0394S', '\u0394H', 'k', '\u0394G_thermo', '\u0394G_k', ]:
        plt.figure(delta)
        delta_list = list()
        for T in T_list:
            delta_list.append(results[T][delta])
        plt.plot(T_list, delta_list, label=delta)
        plt.xlabel('T')
        plt.ylabel(delta.replace('delta', '\u0394'))
        plt.tight_layout()
        if save:
            plt.savefig(save + '_' + delta + '.png')
        plt.show()
    print('\n"THAT IS NOT PORN, THAT IS ART."\n        -- Valeria Jose Boide Trujillo\n')
    return


def visual():
    """>  This function shows on screen all the values calculated and stored in 'data' or 'results'

    Parameters
    ----------
        None

    Returns
    -------
    	None

    """
    df = pd.DataFrame()
    comp = list(data.keys())
    comp.remove('U0')
    T_list = list(data[comp[0]].keys())
    T_list.remove('s_c')
    T_list.remove('param')
    col = ['T', 's_c'] + list(data[comp[0]][T_list[0]].keys()) + ['\u0394U', '\u0394S', '\u0394H', 'k',
                                                                  '\u0394G_thermo', '\u0394G_k', ]
    for T in T_list:
        for compound in comp:
            data[compound][T]['T'] = T
            data[compound][T]['s_c'] = data[compound]['s_c']
            df = pd.concat([df, pd.DataFrame(data[compound][T], index=[compound], columns=col)])
        if T in results:
            df = pd.concat([df, pd.DataFrame(results[T], index=['\u0394'], columns=col)])
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(df)
    return


def save(filename='stat_thermo_results'):
    """>  This function save all the values calculated in a csv and Excel file

    Parameters
    ----------
    filename, optional
        The name of the file

    Returns
    -------
    	None

    """
    df = pd.DataFrame()
    comp = list(data.keys())
    comp.remove('U0')
    T_list = list(data[comp[0]].keys())
    T_list.remove('s_c')
    T_list.remove('param')
    if not T_list:
        print("Data not found")
    col = ['T', 's_c'] + list(data[comp[0]][T_list[0]].keys()) + ['\u0394U', '\u0394S', '\u0394H', 'k',
                                                                  '\u0394G_thermo', '\u0394G_k', ]
    for T in T_list:
        for compound in comp:
            try:
                data[compound][T]['T'] = T
                data[compound][T]['s_c'] = data[compound]['s_c']
                df = pd.concat([df, pd.DataFrame(data[compound][T], index=[compound], columns=col)])
            except:
                print(f"Data not found for {compound} at temperature {T}")
        if T in results:
            df = pd.concat([df, pd.DataFrame(results[T], index=['\u0394'], columns=col)])
    df = df.round(6)
    df.to_csv(filename + '.csv')
    df.to_excel(filename + '.xlsx')
    return 'File saved as ' + filename + ".csv and" + filename + ".xlsx"

