for filename in os.listdir('./local_proteins/'):

    filename = './local_proteins/' + filename
    picklefile = open(filename, 'rb')
    protein = pickle.load(picklefile)
    protein.find_aromatic_rings()

    for outlier_atom in protein.outliers_list:
        