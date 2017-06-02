import numpy as np


def usr_qa_setup(material, log, usr_input):
    '''
    A list of questions to be answered to initialize the CyGutz job.
    '''
    import sys
    print '\n' + " User inputs to initialize the ga job."
    # cut-off distance to determine the rotation group.
    dist_cut = -1.0
    if '-c' in sys.argv:
        dist_cut = float(sys.argv[sys.argv.index('-c') + 1])
        print >> log, " The uniform dist_cut for extracting a centered" + \
            " cluster for symmetry evaluation = ", dist_cut, '\n'
    material.set_sym_dist_cut(dist_cut)

    # Spin symmetry breaking
    spin_polarization = get_usr_input(
        "\n Do you want to BREAK SPIN-SYMMETRY?", ['y', 'n'])
    print >> usr_input, spin_polarization
    material.set_spin_polarization(spin_polarization)

    # Orbital symmetry breaking
    orbital_polarization = get_usr_input(
        "\n Do you want to COMPLETELY BREAK ORBITAL-SYMMETRY?", ['y', 'n'])
    print >> usr_input, orbital_polarization
    material.set_orbital_polarization(orbital_polarization)

    # ASKING USER LIST OF CORRELATED ATOMS AND RELATIVE INFORMATION:
    corr_list = []  # indices correlated atoms
    df_list = []  # information on whether correlated atoms are d or f
    # information on whether correlated atoms have negligible SOC or not
    soc_list = []
    # information on whether correlated atoms have negligible crystal fields
    # or not
    cf_list = []
    if 'y' in spin_polarization:
        Mvec_list = []  # Spin direction list
    else:
        Mvec_list = None
    idx_equivalent_atoms = material.get_myEquivalentAtoms()['equivalent_atoms']
    yn = get_usr_input("\n Symmetrically-equivalent atom indices: " + ''.join(
            "%2d " % (i) for i in idx_equivalent_atoms) +
            "\n (note: '0 0 0 1 1' means 1-3 and 4-5 are two" +
            " inequivalent atoms). \n Accept?", ['y', 'n'])
    print >> usr_input, yn
    if yn == 'n':
        while True:
            string_idx_equivalent_atoms = raw_input(
                " Enter user-defined equivalent atom indices: ")
            print >> usr_input, string_idx_equivalent_atoms
            yn1 = get_usr_input(
                    "\n User-defined equivalent atom indices: " +
                    string_idx_equivalent_atoms + ". Accept?", ['y', 'n'])
            print >> usr_input, yn1
            if yn1 == 'y':
                idx_equivalent_atoms = [
                    int(s) for s in string_idx_equivalent_atoms.split()]
                break
    material.set_idx_equivalent_atoms(idx_equivalent_atoms)

    for i, s in enumerate(material.symbols):
        if i > 0 and idx_equivalent_atoms[i - 1] == idx_equivalent_atoms[i]:
            if 'y' in spin_polarization:
                Mvec_list.append(Mvec_list[-1])
            if 'y' in correlated:
                corr_list.append(i)
                df_list.append(df_list[-1])
                soc_list.append(soc_list[-1])
                cf_list.append(cf_list[-1])
            continue
        print "\n" + ' ' + '-'*12 + "\n atom ", i, " ", s
        correlated = get_usr_input("\n Is this atom correlated?", ['y', 'n'])
        print >> usr_input, correlated
        if correlated == 'y':
            corr_list.append(i)

            df = get_usr_input_combo(
                    "\n Enter correlated shells?", ['s', 'p', 'd', 'f'])
            print >> usr_input, df
            df_list.append(df)

            soc = get_usr_input(
                    "\n Do you want to take into account the SPIN-ORBIT" +
                    " interaction?", ['y', 'n'])
            print >> usr_input, soc
            soc_list.append(soc)

            if 'n' in orbital_polarization:
                cf = get_usr_input(
                    "\n Do you want to take into account the CRYSTAL FIELD" +
                    " effect?", ['y', 'n'])
                print >> usr_input, cf
                cf_list.append(cf)
            else:
                cf_list.append('y')

        if 'y' in spin_polarization:
            if 'n' in correlated:
                Mvec_list.append([0, 0, 0])
            elif 'n' in soc or 'y' in orbital_polarization:
                Mvec_list.append([0, 0, 1])  # z-direction by default.
            else:
                while True:
                    answer = raw_input(
                            "\n Please enter the local spin magnitization" +
                            " direction \n in global coordinate system" +
                            " (e.g., 0 0 1):  ")
                    try:
                        Mvec = [float(v) for v in answer.split()[:3]]
                        print >> usr_input, answer
                        break
                    except:
                        pass
                Mvec = np.array(Mvec)
                Mvec = Mvec / np.linalg.norm(Mvec)
                Mvec_list.append(Mvec)

    material.set_myCorrAtoms(corr_list)
    material.set_myAM(df_list)
    material.set_mySOC(soc_list)
    material.set_myCF(cf_list)
    material.set_Mvec(Mvec_list)


def get_usr_input(message, accept_list):
    while True:
        answer = raw_input(
                message +
                " \n Pick one from [" +
                ', '.join(item for item in accept_list) + "]...")
        if answer not in accept_list:
            print " Please pick an answer in the list!" + \
                    " Make your choice again."
        else:
            break
    return answer


def get_usr_input_combo(message, accept_list):
    while True:
        answer = raw_input(
                message +
                " \n Pick one or combinations separated by blank space" +
                " \n from [" + ', '.join(item for item in accept_list) + "]...")
        if answer_valid(answer, accept_list):
            break
    return answer


def answer_valid(answer, accept_list):
    answer_list = answer.split()
    for ans in answer_list:
        if ans not in accept_list:
            return False
    return True


def inp_ga_init_dmft(soc_list, log):
    '''
    Generate the inputs for running ga_init_dmft.py.
    '''
    import sys
    if 'vasp' in sys.argv:
        return
    print " Please run ga_init_dmft.py with parameters given in", log.name
    print >> log, '\n' + \
        " Please run ga_init_dmft.py with the following parameters."
    if 'y' in soc_list:
        qsplit = 4
    else:
        qsplit = 3
    print >> log, '\n' + \
        ' QUESTION: Specify qsplit for each correlated orbital (default = 0)'
    print >> log, ' ANSWER: ', qsplit
    print >> log, '\n' + \
        ' QUESTION: Do you want to group any of these orbitals' + \
        'into cluster-DMFT problems? (y/n)'
    print >> log, ' ANSWER: ', 'n'
    answer = ' '.join([str(i + 1) for i in range(len(soc_list))])
    print >> log, '\n' + \
        ' QUESTION: Enter the correlated problems' + \
        'forming each unique correlated problem,' + \
        'separated by spaces (ex: 1,3 2,4 5-8)'
    print >> log, ' ANSWER: ' + answer
    print >> log, '\n' + ' QUESTION: Broken symmetry run? (y/n)'
    print >> log, ' ANSWER:', 'n'
    print >> log, '\n' + ' QUESTION: Is this a spin-orbit run? (y/n)'
    if 'y' in soc_list:
        print >> log, ' ANSWER:', 'y'
    else:
        print >> log, ' ANSWER:', 'n'
    print >> log, '\n' + '\n'
