def read_line_catalog(input_catalog):

    """
    (xpos_atoms, labels_atoms)  = read_line_catalog('auxilarary_data/line_catalogs/edibles_linelist_atoms.csv')
    """

    x = np.genfromtxt(input_catalog, dtype=None, skip_header=1, delimiter=',', unpack=True)
    xpos  = x['f1']
    labels= x['f2']


    """ alternative:
    xpos,labels = np.genfromtxt(input_catalog, usecols=(0,1), skip_header=1, dtype=None, delimiter=',')
    """
    
    return xpos, labels
