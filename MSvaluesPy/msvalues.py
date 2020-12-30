import pkg_resources
import numpy as np, pandas as pd
import os, sys, scipy.linalg.blas, scipy.stats
from cvxopt import matrix, solvers
import matplotlib.pyplot as plt
from collections import Counter


# sets up family-specific marker effects for Mendelian sampling
def makemems(geno, me):
    matrix_ms = geno
    matrix_ms.astype(np.float64)
    matrix_ms = np.where(matrix_ms < 3., 0., matrix_ms)
    for iiii in range(matrix_ms.shape[1]):
        matrix_ms[:,iiii] = np.where(matrix_ms[:,iiii] == 3., (me[iiii]), matrix_ms[:,iiii])
        matrix_ms[:,iiii] = np.where(matrix_ms[:,iiii] == 4., (me[iiii]*-1), matrix_ms[:,iiii])
    return matrix_ms


# converts covariance to correlation matrix
def correlation_from_covariance(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation

# Import example input file 1: marker map and trait effects for a single trait
def testinputfile1st():
    cwd = os.getcwd()
    stream = pkg_resources.resource_stream(__name__, 'example_data/inputfile1ST.txt')
    data = pd.read_csv(stream, sep = " ")
    os.chdir(cwd)
    return data

# Import example input file 2: ID, sex, phased genotype data coded 1, 2, 3, 4 (for single sex)
def testinputfile2ss():
    cwd = os.getcwd()
    stream = pkg_resources.resource_stream(__name__, 'example_data/inputfile2SS.txt')
    data = np.loadtxt(stream, dtype=str)
    os.chdir(cwd)
    return data

# Import example input file 1: marker map and trait effects for a multiple trait
def testinputfile1mt():
    cwd = os.getcwd()
    stream = pkg_resources.resource_stream(__name__, 'example_data/inputfile1MT.txt')
    data = pd.read_csv(stream, sep = " ")
    os.chdir(cwd)
    return data

# Import example input file 2: ID, sex, phased genotype data coded 1, 2, 3, 4 (for multiple sex)
def testinputfile2ms():
    cwd = os.getcwd()
    stream = pkg_resources.resource_stream(__name__, 'example_data/inputfile2MS.txt')
    data = np.loadtxt(stream, dtype=str)
    os.chdir(cwd)
    return data

# checks the input data files to correctness
def fileschecks(inputfile1, inputfile2, no_traits, index_wt):
    inputfile1 = inputfile1
    inputfile2 = inputfile2
    no_traits = no_traits #no_traits
    index_wt = np.array(index_wt) #inde weights
    # check: ensures number of traits is an integer
    if (type(no_traits) != int):
        sys.exit('Number of traits must be an integer')
    # check: ensures number of traits match size of index weights
    if (no_traits != index_wt.size):
        sys.exit('Index weights do not match no of trait(s)')
    no_markers = inputfile1.shape[0] # no of markers = length of file
    no_chromosomes = max(inputfile1.iloc[:,1]) # no of chromosomes
    trait_names = inputfile1.columns[3:len(inputfile1.columns)] # traits names
    # check: ensures number of traits match length of trait names/marker effects in inputfile 1
    if (no_traits != trait_names.size):
        sys.exit("No. of traits does not match data in marker map and marker effects dataframe")
    #check: to ensure the first column is appropriately numbered. It is used by the progra for subsetting
    compare = np.array(inputfile1.iloc[:,0]) == np.sort(np.array(inputfile1.iloc[:,0]))
    if (compare.any() == False):
        sys.exit("The first column of marker map and marker effects dataframe should be numbered from 1 to the total number of markers")
    del compare
    # ensures the marker position is ordered
    for chr in range(no_chromosomes):
        mp = np.array(inputfile1.iloc[:, 2][inputfile1.iloc[:, 1] == (1+chr)])
        smp = np.sort(sorted(mp))
        compare = mp == smp
        if (compare.any() == False):
            sys.exit("Marker map is faulty on chromosome", (1+chr))

    no_individuals = inputfile2.shape[0]
    ID = inputfile2[:,0] # store animal ID
    ID = ID.astype(str) # converts the IDs to strings
    no_individuals = inputfile2.shape[0] # number of individuals (no. of rows of input file 2)
    # check for duplicated ID of animals
    if (len(set(ID)) != no_individuals):
            sys.exit("Duplicates exist in the IDs of individuals (first column of id, sex and genotypes array)")
    sex = inputfile2[:,1] # store animal ID
    sex = sex.astype(int) # converts the sex to intger
    #check: to ensure the sex data is coded 1, 2
    if (min(sex) < 1 or max(sex) > 2):
        sys.exit("Check values for sex. Sex should be coded using 1 for male and 2 for female")

    inputfile2 = inputfile2[:,2:inputfile2.shape[1]]  # store genotype data
    inputfile2 = pd.DataFrame(inputfile2)
    inputfile2 = inputfile2.values.astype(int) # converts to int
    if (inputfile2.max() > 4):
        sys.exit("Genotypes in id, sex and genotypes array contain value(s) greater than 4")
    elif (inputfile2.min() < 1):
        sys.exit("Genotypes in id, sex and genotypes array contain value(s) less than 1")
    if (no_markers != inputfile2.shape[1]):
        sys.exit("Discrepancy between no. markers in (marker map and effects dataframe) and 2 (id, sex and genotypes array)")
    print('Input files or data passed the file check test!')
    print("Number of individuals: ",no_individuals)
    print("Number of chromosomes: ", no_chromosomes)
    print("Total no. markers: ", no_markers)
    print("Number of trait(s): ", no_traits)
    print("Trait name(s) and Index weight(s)")
    if (no_traits == 1):
        for i in range(no_traits):
            print(trait_names[i], ": ", index_wt)
    elif (no_traits > 1):
        for i in range(no_traits):
            print(trait_names[i], ": ", index_wt[i])


# sets up chromosome-wise population covariance matrix
def elem_cor(list, mp, chr, mposunit):
    if (mposunit == "bp" or mposunit == "Bp" or mposunit == "bP" or mposunit == "BP"):
        mp = mp/1000000
    tmp = np.exp(-2*(np.abs(mp - mp[:, None])*0.01))/4
    list[0].append(tmp)
    return list


# sets up population covariance matrix as a list of chromosome-wise population covariance matrices
def popcovmat(inputfile1, mposunit):
    mposunit = mposunit
    inputfile1 = inputfile1
    no_chromosomes = max(inputfile1.iloc[:,1]) # no of chromosomes
    list = [] #list stores chromosome-wise covariance matrix
    list.append([])
    # print("Chromosomes processed:")
    for chr in range(no_chromosomes):
        # select marker positions (mp) correcponding to each chromosome
        mp = np.array(inputfile1.iloc[:, 2][inputfile1.iloc[:, 1] == (1+chr)])
        elem_cor(list, mp, chr, mposunit)
    return list

# prints the progress of a function
def progr(no_individuals, xx):
    prog = np.ceil(np.quantile(np.arange(1, no_individuals, 1).tolist(), [.25, .5, .75, 1]))
    prog = prog.astype(int)
    if (xx == 0):
        print("Progress (%): ", end = "", flush=True)
    elif (xx == prog[0]):
        print("25", " ", end = "", flush=True)
    elif (xx == prog[1]):
        print("50", " ", end = "", flush=True)
    elif (xx == prog[2]):
        print("75", " ", end = "", flush=True)
    elif (xx == prog[3]):
        print("100", " ", flush=True)

# list to store trait-specific matrices
def traitspecmatrices(inputfile1, inputfile2, no_traits):
    trait_effects = inputfile1.iloc[:, 3:inputfile1.shape[1]] # traits effects
    trait_effects = trait_effects.to_numpy() # convert dataframe to array    
    genomat = inputfile2[:,2:inputfile2.shape[1]]  # store genotype data
    genomat = pd.DataFrame(genomat)
    genomat = genomat.values.astype(int) # converts to int
    genomat = np.array(genomat)
    psvlist = [] #list stores trait-specific matrices
    psvlist.append([])
    for i in range(no_traits):
        matrix_ms = makemems(genomat, trait_effects[:,i])
        psvlist[0].append(matrix_ms)
    return psvlist

# obtains column names for MSVMSV dataframe
def colnamesfordf(no_traits, trait_names):
    tn = np.zeros((no_traits, no_traits), 'U20')
    tn = np.chararray(tn.shape, itemsize=20)
    for i in range(no_traits):
        for t in range(no_traits):
            if (i == t):
                tn[i,t] = str(trait_names[i])
            elif (i != t):
                tn[i,t] = "{}_{}".format(trait_names[i],trait_names[t])
    colnamesms = tn[np.tril_indices(no_traits)]
    return colnamesms

# MRM' using np.matmul
def MRMmult(temp, CovMat):
    return np.matmul(np.matmul(temp, CovMat), temp.transpose())        

# # derives Mendelian sampling co(variance) and aggregate breeding value if multiple traits
def MSvarcov(inputfile1, inputfile2, CovMat, index_wt, progress = True):
    ID = inputfile2[:,0] # store animal ID
    ID = ID.astype(str) # converts the IDs to strings
    sex = inputfile2[:,1] # store animal ID
    sex = sex.astype(int) # converts the IDs to strings
    no_individuals = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    no_chromosomes = max(inputfile1.iloc[:,1]) # no of chromosomes
    trait_names = inputfile1.columns[3:len(inputfile1.columns)] # traits names
    no_traits = trait_names.size
    
    psvlist = traitspecmatrices(inputfile1, inputfile2, no_traits) # list to store trait-specific matrices    

    # derives column names for dataframe for storing results
    colnamesms = colnamesfordf(no_traits, trait_names)
    
    # dataframe to save Mendelian sampling (co)variance and aggregate breeding value if multu
    MSVMSC = np.zeros((no_individuals, no_traits*(no_traits-1)+1))
    for i in range(no_individuals): #oop over no of individuals
        MSC = np.zeros((no_traits, no_traits))
        for chr in range(no_chromosomes): #Loop over the number of chromosomes
            snpindex = np.array(inputfile1.iloc[:, 0][inputfile1.iloc[:, 1] == (1+chr)])-1 #marker index
            temp = np.zeros((no_traits, len(snpindex)))      
            for tr in range(no_traits):
                temp[tr, :] = psvlist[0][tr][i,snpindex]
            MSC = MSC + MRMmult(temp, CovMat[0][chr])
        if (no_traits == 1):
            MSVMSC[i,0] = MSC
        elif (no_traits > 1):
            iu1 = np.tril_indices(no_traits) #extracts only low triangular matrix
            MSVMSC[i, 0:(no_traits*(no_traits-1))] = MSC[iu1] #flatten the matrix in predetermined order
            MSVMSC[i, (no_traits*(no_traits-1))] = MRMmult(index_wt, MSC)

        if progress:
            # prints progress in percentage
            progr(no_individuals, i)     
    MSVMSC=pd.DataFrame(MSVMSC) #convert matrix to dataframe
    if (no_traits == 1):
        MSVMSC.columns = trait_names
    elif (no_traits > 1):
        colnamesms = np.concatenate((colnamesms, "Agg BV"), axis=None)
        MSVMSC.columns = colnamesms
    MSVMSC.insert(0, "Animal_ID", ID, True) # insert animal ID in the first col of the dataframe
    MSVMSC.insert(1, "Sex", sex, True) # insert animal ID in the first col of the dataframe
    return MSVMSC



# sets up family-specific marker effects for GEBV estimation
def makemebv(inputfile2, me):
    matrix_me = inputfile2[:,2:inputfile2.shape[1]]  # store genotype data
    matrix_me = matrix_me.astype(np.float64)
    matrix_me = np.where(matrix_me > 2., 0., matrix_me)
    for iiii in range(matrix_me.shape[1]):
        matrix_me[:,iiii] = np.where(matrix_me[:,iiii] == 1., (me[iiii]), matrix_me[:,iiii])
        matrix_me[:,iiii] = np.where(matrix_me[:,iiii] == 2., (me[iiii]*-1), matrix_me[:,iiii])
    return matrix_me


# calculates GEBV
def CalcGEBV(inputfile1, inputfile2, index_wt):
    ID = inputfile2[:,0] # store animal ID
    ID = ID.astype(str) # converts the IDs to strings
    sex = inputfile2[:,1] # store animal ID
    sex = sex.astype(int) # converts the IDs to strings
    no_individuals = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    trait_names = inputfile1.columns[3:len(inputfile1.columns)] # traits names
    trait_effects = inputfile1.iloc[:, 3:inputfile1.shape[1]] # traits effects
    trait_effects = trait_effects.to_numpy() # convert dataframe to array
    no_traits = trait_names.size

    # calculates GEBV
    if (no_traits == 1):
        GEBV = np.zeros((no_individuals, no_traits))
        matrix_me = makemebv(inputfile2, trait_effects[:,0])
        GEBV[:, 0] = matrix_me.sum(axis=1)
        GEBV=pd.DataFrame(GEBV)
        GEBV.columns = trait_names
    elif (no_traits > 1):
        GEBV = np.zeros((no_individuals, no_traits+1))
        for i in range(no_traits):
            matrix_me = makemebv(inputfile2, trait_effects[:,i])
            GEBV[:, i] = matrix_me.sum(axis=1)
            GEBV[:, no_traits] = GEBV[:, no_traits] + index_wt[i]*GEBV[:, i]
        GEBV=pd.DataFrame(GEBV) #convert matrix to dataframe
        colnames = np.concatenate((trait_names, "Overall"), axis=None)
        GEBV.columns = colnames
    GEBV.insert(0, "Animal_ID", ID, True) # insert animal ID in the first col of the dataframe
    GEBV.insert(1, "Sex", sex, True) # insert animal ID in the second col of the dataframe
    return GEBV

# calculates PBTI using a threshold top x% of GEBV
def CalcProb(inputfile1, inputfile2, MSVMSC, index_wt, thresh):
    ID = inputfile2[:,0] # store animal ID
    ID = ID.astype(str) # converts the IDs to strings
    sex = inputfile2[:,1] # store animal ID
    sex = sex.astype(int) # converts the IDs to strings
    no_individuals = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    trait_names = np.array(inputfile1.columns[3:len(inputfile1.columns)]) # traits names
    trait_effects = inputfile1.iloc[:, 3:inputfile1.shape[1]] # traits effects
    trait_effects = trait_effects.to_numpy() # convert dataframe to array
    no_traits = trait_names.size

    GEBV = CalcGEBV(inputfile1, inputfile2, index_wt)
    #Prob to breed top individuals
    if (no_traits == 1):
        Probdf = np.zeros((no_individuals, no_traits))
        t = np.quantile(GEBV.iloc[:,(0+2)], q = 1-thresh)
        Probdf[:,0] = 1 - scipy.stats.norm.cdf(t, loc = GEBV.iloc[:,(0+2)], scale = np.sqrt(MSVMSC.iloc[:,0+2]))
        Probdf=pd.DataFrame(Probdf)
        Probdf.columns = trait_names
    elif (no_traits > 1):
        # multi-trait prob
        colnamesms = colnamesfordf(no_traits, trait_names)        
        colnamesms = np.concatenate((colnamesms, "Agg BV"), axis=None)
        t = np.quantile(GEBV.iloc[:,(no_traits+2)], q = 1-thresh)
        Probdf = np.zeros((no_individuals, no_traits+1))
        traitindex = np.arange(colnamesms.shape[0])[np.in1d(colnamesms, trait_names)]
        for i in range(no_traits):
            t = np.quantile(GEBV.iloc[:,(i+2)], q = 1-thresh)
            Probdf[:,i] = 1 - scipy.stats.norm.cdf(t, loc = GEBV.iloc[:,(i+2)], scale = np.sqrt(MSVMSC.iloc[:,(traitindex[i])+2]))
            # multi-trait prob
        t = np.quantile(GEBV.iloc[:,(no_traits+2)], q = 1-thresh)
        Probdf[:,no_traits] = 1 - scipy.stats.norm.cdf(t, loc = GEBV.iloc[:,(no_traits+2)], scale = np.sqrt(MSVMSC.iloc[:,(no_traits*(no_traits-1))+2]))
        Probdf=pd.DataFrame(Probdf) #convert matrix to dataframe
        colnames = np.concatenate((trait_names, "Overall"), axis=None)
        Probdf.columns = colnames
    Probdf.insert(0, "Animal_ID", ID, True) # insert animal ID in the first col of the dataframe
    Probdf.insert(1, "Sex", sex, True) # insert animal ID in the second col of the dataframe
    return Probdf

# calculates RPTA if selection intensity is known
def CalcRPTA(inputfile1, inputfile2, MSVMSC, index_wt, selint): 
    ID = inputfile2[:,0] # store animal ID
    ID = ID.astype(str) # converts the IDs to strings
    sex = inputfile2[:,1] # store animal ID
    sex = sex.astype(int) # converts the IDs to strings
    no_individuals = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    trait_names = np.array(inputfile1.columns[3:len(inputfile1.columns)]) # traits names
    trait_effects = inputfile1.iloc[:, 3:inputfile1.shape[1]] # traits effects
    trait_effects = trait_effects.to_numpy() # convert dataframe to array
    no_traits = trait_names.size
    GEBV = CalcGEBV(inputfile1, inputfile2, index_wt)
    #RPTA

    if (no_traits == 1):
        RPTAdf = np.zeros((no_individuals, no_traits))
        RPTAdf[:,0] = (GEBV.iloc[:,(0+2)]/2) + np.sqrt(MSVMSC.iloc[:,0+2])*selint
        RPTAdf=pd.DataFrame(RPTAdf)
        RPTAdf.columns = trait_names
    elif (no_traits > 1):
        # multi-trait prob
        colnamesms = colnamesfordf(no_traits, trait_names)
        colnamesms = np.concatenate((colnamesms, "Agg BV"), axis=None)

        RPTAdf = np.zeros((no_individuals, no_traits+1))
        traitindex = np.arange(colnamesms.shape[0])[np.in1d(colnamesms, trait_names)]
        for i in range(no_traits):
            RPTAdf[:,i] = (GEBV.iloc[:,(i+2)]/2) + np.sqrt(MSVMSC.iloc[:,(traitindex[i]+2)])*selint
            # multi-trait RPTA
        RPTAdf[:,no_traits] = (GEBV.iloc[:,(no_traits+2)]/2) + np.sqrt(MSVMSC.iloc[:,(no_traits*(no_traits-1))+2])*selint
        RPTAdf=pd.DataFrame(RPTAdf)
        colnames = np.concatenate((trait_names, "Overall"), axis=None)
        RPTAdf.columns = colnames
    RPTAdf.insert(0, "Animal_ID", ID, True) # insert animal ID in the first col of the dataframe
    RPTAdf.insert(1, "Sex", sex, True) # insert animal ID in the second col of the dataframe
    return RPTAdf

# Calculates selection criteria (GEBV, PBTI, or RPTA) using gametic approach
def SelStratGam(selstrat, inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh):
    if (selstrat == 'GEBV'):
         data = CalcGEBV(inputfile1, inputfile2, index_wt)
    elif (selstrat == 'PBTI'):
        data = CalcProb(inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh)
    elif (selstrat == 'RPTA'):
        data = CalcRPTA(inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh)
    return data

# creates matrix for additive effects of aggregate genotype
def agggenmse(no_individuals, no_markers, Mlist, index_wt):
    MMfinal =  np.empty((no_individuals, no_markers))
    for i in range(no_individuals):
        tempmat1 = np.zeros((index_wt.size, no_markers))
        for tr in range(index_wt.size):
            tempmat1[tr,:] = Mlist[0][tr][i,:]
        MMfinal[i,:] = np.matmul(index_wt.transpose(), tempmat1)
    return MMfinal

# MRM' using scipy.linalg.blas.dgemm
def dgemmMRM(temp, CovMat):
    temp1111 = scipy.linalg.blas.dgemm(alpha = 1.0, a = temp, b = CovMat) #M*R
    covtemp = abs(scipy.linalg.blas.dgemm(alpha = 1.0, a = temp1111, b = temp.transpose()))  #((M*R)*M'))
    return covtemp
    
# Derives similarity matrices using gametic approach    
def SimMatGam(inputfile1, inputfile2, index_wt, CovMat, chrinterest, stdsim = True, progress=True):
    chrinterest = chrinterest
    if ('all' in chrinterest):
        chrinterest = 'all'
    elif ('none' in chrinterest):
        chrinterest = 'none'
    else:
        chrinterest = [int(i) for i in chrinterest]
        chrinterest = np.array(chrinterest)

    no_individuals = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    no_markers = inputfile1.shape[0] # no of markers = length of file
    no_chromosomes = max(inputfile1.iloc[:,1]) # no of chromosomes
    trait_names = inputfile1.columns[3:len(inputfile1.columns)] # traits names
    no_traits = trait_names.size
    
    psvlist = traitspecmatrices(inputfile1, inputfile2, no_traits) # list to store trait-specific matrices    

    # The following estimates similarity matrices
    #loop over each trait
    for i in range(no_traits):
        # create a matrix of zeros with rows and cols = no of inds
        cov_bet_indxxxxxxxx = np.zeros((no_individuals, no_individuals)) #stores trait-specific covariance btw inds
        for chr in range(no_chromosomes): #Loop over the number of chromosomes
            if progress:
                if (chr == 0):
                    print("Processing  ", trait_names[i])
            snpindex = np.array(inputfile1.iloc[:, 0][inputfile1.iloc[:, 1] == (1+chr)])-1 #marker index
            #matrix multiplication using BLAS
            covtempxxxxxxxxxxxxxxxxx = abs(dgemmMRM(psvlist[0][i][:,snpindex], CovMat[0][chr]))
            cov_bet_indxxxxxxxx = cov_bet_indxxxxxxxx + covtempxxxxxxxxxxxxxxxxx #sums up chr specific covariances
            if (type(chrinterest) == str):
                if (chrinterest == 'all'):
                    chrfile = "{}/Similarity matrix btw ind_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                    np.save(chrfile, covtempxxxxxxxxxxxxxxxxx) # writes sim matrices for chromosomes of interest to file
            elif (chr+1 in chrinterest):
                chrfile = "{}/Similarity matrix btw ind_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                np.save(chrfile, covtempxxxxxxxxxxxxxxxxx) # writes sim matrices for chromosomes of interest to file
            if stdsim:
                if (type(chrinterest) == str):
                    if (chrinterest == 'all'):
                        corrtempxxxxxxxxxxx = correlation_from_covariance(covtempxxxxxxxxxxxxxxxxx)
                        chrfilec = "{}/Standardized similarity matrix btw ind_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                        np.save(chrfilec, corrtempxxxxxxxxxxx) # writes std sim matrices for chromosomes of interest to file
                        del corrtempxxxxxxxxxxx
                elif (chr+1 in chrinterest):
                    corrtempxxxxxxxxxxx = correlation_from_covariance(covtempxxxxxxxxxxxxxxxxx)
                    chrfilec = "{}/Standardized similarity matrix btw ind_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                    np.save(chrfilec, corrtempxxxxxxxxxxx) # writes std sim matrices for chromosomes of interest to file
                    del corrtempxxxxxxxxxxx
            del covtempxxxxxxxxxxxxxxxxx
            #outputs progress in terms of trait name and chr processed
            if progress:
                progr(no_chromosomes, chr)

    if (no_traits == 1):
        mat = cov_bet_indxxxxxxxx
        if stdsim:
            mat = correlation_from_covariance(mat)
    elif (no_traits > 1):
        if progress:
            print('Creating similarity matrix based on aggregate genotype')
        tempmat1 = agggenmse(no_individuals, no_markers, psvlist, index_wt)
        mat = np.zeros((no_individuals, no_individuals)) # stores overall covariance btw inds.
        # loop over chromososomes
        for chr in range(no_chromosomes):
            snpindex = np.array(inputfile1.iloc[:, 0][inputfile1.iloc[:, 1] == (1+chr)])-1
            covtempxxxxxxxxxxxxxxxxx = abs(dgemmMRM(tempmat1[:,snpindex], CovMat[0][chr])) #((M*R)*M')
            mat = mat + covtempxxxxxxxxxxxxxxxxx
            del covtempxxxxxxxxxxxxxxxxx

            if progress:
                progr(no_chromosomes, chr)

        if stdsim:
            mat = correlation_from_covariance(mat)
    return mat


# Plots and saves EF
def PlotEF(opt_portfolios, Similaritymat, returns, no_individuals):    
    minvar = opt_portfolios[:,0]
    index_min = np.argmin(minvar)
    minsd_er = opt_portfolios[index_min,1]
    efficient = opt_portfolios[:,1] >= minsd_er
    aa = np.where(efficient == 0)
    bb = np.where(efficient == 1)
    Allind = np.zeros((no_individuals, 2))
    Allind[:,0], Allind[:,1] = np.diag(Similaritymat), returns
    # Plot EF
    plt.figure(1)
    plt.scatter(Allind[:,0], Allind[:,1], color = 'grey', marker = ".", s = 8)
    plt.scatter(opt_portfolios[aa,0], opt_portfolios[aa,1], color = 'red', marker = "o", s = 8)
    plt.scatter(opt_portfolios[bb,0], opt_portfolios[bb,1], color = 'blue',  marker = "o", s = 8)
    plt.scatter(opt_portfolios[index_min,0], opt_portfolios[index_min,1], color = 'blue', s = 9,  marker = "s")
    plt.scatter(opt_portfolios[np.argmin(opt_portfolios[:,0]/opt_portfolios[:,1]),0],
            opt_portfolios[np.argmin(opt_portfolios[:,0]/opt_portfolios[:,1]),1], color = 'blue', s = 10,  marker = "D")
    plt.xlabel("Variance")
    plt.ylabel("Expected return")
    plt.title("Optimum mate allocation using efficient frontier")
    #plt.show()
    plt.savefig('Efficient Frontier_gametic approach.pdf')
    
# Optimizes mate allocation   
def OptiMateAllocGam(data, Similaritymat, maxalloc, progress = True):
    no_individuals = data.shape[0]
    result = np.max(data.iloc[:,1]) != np.min(data.iloc[:,1])
    if result:
        sys.exit('Individuals are not of the same sex. Use zygotic approach')

    # Expected return on progeny of an individual
    if (data.shape[1] == 3):
        returns = data.iloc[:,2]
    elif (data.shape[1] > 3):
        returns = data['Overall']

    # Optimizing mate allocation
    muhat = np.linspace(start=min(returns), stop=max(returns), num=50)
    sln = np.zeros((len(muhat), no_individuals))
    opt_portfolios = np.zeros((len(muhat), 2))
    Similaritymat = matrix(Similaritymat)  # Covariance matrix
    q = matrix(np.zeros((no_individuals, 1)))
    r_avg = np.zeros((1, no_individuals))
    r_avg[0,:] = returns   #selstrat values

    # Capturing inequality constraints Gx >= h
    # captures the constraints (r_avg'x >= muhat[i]) and (x >= lb and x <= ub)
    G = matrix(np.vstack((-r_avg, -1.0*np.eye(no_individuals) , np.eye(no_individuals))))
    L=np.array([-0.0]*no_individuals) #lower bound of 0.0
    # upper bound. to be provided from param file
    U=np.array([maxalloc]*no_individuals)
    # Capturing equality constraint (Ax == b)
    A = matrix(1.0, (1,no_individuals))
    b = matrix(1.0)
    
    if progress:
        print('Optimizing portfolios')
    for i in range(len(muhat)):
        try:
            h=matrix(np.transpose(np.hstack((-muhat[i], L, U)).reshape((1,1+(2*no_individuals)))))
            solvers.options['show_progress'] = False
            sol = solvers.qp(Similaritymat, q, G, h, A, b)
            sln[i,:] = np.array(sol['x']).transpose()
            opt_portfolios[i,0] = np.matmul(np.matmul(sln[i,:].transpose(), Similaritymat), sln[i,:])
            opt_portfolios[i,1] = np.dot(sln[i,:].transpose(), returns)

            if progress:
                progr(len(muhat), i)
        except Exception as err:
            exception_type = type(err).__name__
    if progress:
        if (maxalloc != 1):
            print("100", " ", flush=True)
        print('Mate allocation optimized')
    sln = sln[~np.all(sln == 0, axis=1)]
    slnfile = "{}/Optimum allocation of weights to individuals.csv".format(os.getcwd()) # output file
    sln = pd.DataFrame(sln)
    sln.columns = data.iloc[:,0]
    sln.to_csv(slnfile, header = True, index=False) # prints portfolio weights
    
    opt_portfolios = opt_portfolios[~np.all(opt_portfolios == 0, axis=1)]
    
    PlotEF(opt_portfolios, Similaritymat, returns, no_individuals)
    
    opt_portfolios = pd.DataFrame(opt_portfolios)
    opt_portfolios.columns = ['Risk', 'Returns']
    optfile = "{}/Expected risk and returns on portfolios.csv".format(os.getcwd()) # output file
    opt_portfolios.to_csv(optfile, header = True, index=False) # prints portfolio risks and returns
    return opt_portfolios, sln


def potentialparents(selstrat, data, no_traits, selintmale, selintfemale, trait_names):
    if (no_traits == 1):
        datamale = data[data.iloc[:,1] == 1]
        pos = data.index.values[data.iloc[:,1] == 1]
        datamale.insert(0, "pos", pos, True)
        no_sire = int(datamale.shape[0] * selintmale) # 10% of the number of sire
        datamale = datamale.sort_values(by=[trait_names], ascending=False).iloc[0:no_sire,:] #
        datafemale = data[data.iloc[:,1] == 2]
        pos = data.index.values[data.iloc[:,1] == 2]
        datafemale.insert(0, "pos", pos, True)
        no_dam = int(datafemale.shape[0] * selintfemale)
        datafemale = datafemale.sort_values(by=[trait_names], ascending=False).iloc[0:no_dam,:]
    elif (no_traits > 1):
        datamale = data[data.iloc[:,1] == 1]
        pos = data.index.values[data.iloc[:,1] == 1]
        datamale.insert(0, "pos", pos, True)
        no_sire = int(datamale.shape[0] * selintmale) # 10% of the number of sire
        datamale = datamale.sort_values(by=['Overall'], ascending=False).iloc[0:no_sire,:] #
        datafemale = data[data.iloc[:,1] == 2]
        pos = data.index.values[data.iloc[:,1] == 2]
        datafemale.insert(0, "pos", pos, True)
        no_dam = int(datafemale.shape[0] * selintfemale) # 50% of the number of sire
        datafemale = datafemale.sort_values(by=['Overall'], ascending=False).iloc[0:no_dam,:] #
        
    sire = datamale.iloc[:,0]
    dam = datafemale.iloc[:,0]
    matlist = np.array(np.meshgrid(sire, dam)).T.reshape(-1,2)
    SireID = datamale.iloc[:,1]
    DamID = datafemale.iloc[:,1]
    IDS = np.array(np.meshgrid(SireID, DamID)).T.reshape(-1,2)
    if (selstrat =='GEBV'):
        matingdata = pd.DataFrame(index=range(matlist.shape[0]),columns=range(6))
    elif (selstrat == 'RPTA' or selstrat == 'PBTI'):
        matingdata = pd.DataFrame(index=range(matlist.shape[0]),columns=range(7))
    matingdata.iloc[:,[0,1]] = IDS
    matingdata.iloc[:,[2,3]] = matlist
    return matingdata


# Calculates selection criteria (GEBV, PBTI, or RPTA) using zygotic approach and saves 
# Mendelian sampling (co)variance for potential parents to file
def SelStratZyg(selstrat, inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh, selintmale, selintfemale):
    trait_names = np.array(inputfile1.columns[3:len(inputfile1.columns)]) # traits names
    no_traits = trait_names.size
    if (selstrat == "GEBV"):
        GEBV = CalcGEBV(inputfile1, inputfile2, index_wt)
        matingdata = potentialparents(selstrat, GEBV, no_traits, selintmale, selintfemale, trait_names)
        if (no_traits == 1):
            matingdata.iloc[:, 4] = np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 5] = np.array(MSVMSC.iloc[:,0][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:,0][matingdata.iloc[:,3]])
        elif (no_traits > 1):
            matingdata.iloc[:, 4] = np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 5] = np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,3]])

        idfxxxxxxxxxxx = np.unique(matingdata.iloc[:,3])
        newmatingmatxxxxxxxxxxx = pd.DataFrame(index=range(len(idfxxxxxxxxxxx)),columns=range(6))
        for mm in range(len(idfxxxxxxxxxxx)):
            axxxxxxxxxxxx = matingdata.loc[matingdata.iloc[:,3] == idfxxxxxxxxxxx[mm]]
            topsire = np.array(axxxxxxxxxxxx.iloc[:,2])
            newmatingmatxxxxxxxxxxx.iloc[mm, :] = axxxxxxxxxxxx.iloc[np.argmax(axxxxxxxxxxxx.iloc[:,4]),:]
            norepssire = Counter(newmatingmatxxxxxxxxxxx.iloc[:,2])
            for nr in range(len(topsire)):
                if (norepssire[topsire[nr]] <= 49):
                    newmatingmatxxxxxxxxxxx.iloc[mm, :] = np.array(axxxxxxxxxxxx[axxxxxxxxxxxx.iloc[:,2]==topsire[nr]])
                    break
        
        GEBVfile = "{}/{}_gametic approach.csv".format(os.getcwd(), selstrat) # output file
        GEBV.to_csv(GEBVfile, header = True, index=False)
        matingdata = newmatingmatxxxxxxxxxxx
        del newmatingmatxxxxxxxxxxx, axxxxxxxxxxxx, idfxxxxxxxxxxx
        matingdata.columns = ['SireID', 'DamID', 'SireIndex', 'DamIndex', 'GEBV', 'MSV']
        matingdatafile = "{}/Mating data (Zygotic).csv".format(os.getcwd()) # output file
        matingdata.to_csv(matingdatafile, header = True, index=False) # prints GEBV of breeding top inds to file


    elif (selstrat == 'PBTI'):
        Probdf = CalcProb(inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh)
        matingdata = potentialparents(selstrat, Probdf, no_traits, selintmale, selintfemale, trait_names)

        GEBV = CalcGEBV(inputfile1, inputfile2, index_wt)
        if (no_traits == 1):
            matingdata.iloc[:, 5] = np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 6] = np.array(MSVMSC.iloc[:,0+2][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:,0+2][matingdata.iloc[:,3]])
            t = np.quantile(matingdata.iloc[:, 5], q = 1-selinorthresh)
            matingdata.iloc[:, 4] = 1 - scipy.stats.norm.cdf(t, loc = matingdata.iloc[:, 5], scale = np.sqrt(matingdata.iloc[:, 6]))
        elif (no_traits > 1):
            matingdata.iloc[:, 5] = np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 6] = np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,3]])
            t = np.quantile(matingdata.iloc[:, 5], q = 1-selinorthresh)
            matingdata.iloc[:, 4] = 1 - scipy.stats.norm.cdf(t, loc = matingdata.iloc[:, 5], scale = np.sqrt(matingdata.iloc[:, 6]))


        idfxxxxxxxxxxx = np.unique(matingdata.iloc[:,3])
        newmatingmatxxxxxxxxxxx = pd.DataFrame(index=range(len(idfxxxxxxxxxxx)),columns=range(7))
        for mm in range(len(idfxxxxxxxxxxx)):
            axxxxxxxxxxxx = matingdata.loc[matingdata.iloc[:,3] == idfxxxxxxxxxxx[mm]]
            topsire = np.array(axxxxxxxxxxxx.iloc[:,2])
            newmatingmatxxxxxxxxxxx.iloc[mm, :] = axxxxxxxxxxxx.iloc[np.argmax(axxxxxxxxxxxx.iloc[:,4]),:]
            norepssire = Counter(newmatingmatxxxxxxxxxxx.iloc[:,2])
            for nr in range(len(topsire)):
                if (norepssire[topsire[nr]] <= 49):
                    newmatingmatxxxxxxxxxxx.iloc[mm, :] = np.array(axxxxxxxxxxxx[axxxxxxxxxxxx.iloc[:,2]==topsire[nr]])
                    break
        
        Probdffile = "{}/{}_gametic approach.csv".format(os.getcwd(), selstrat) # output file
        Probdf.to_csv(Probdffile, header = True, index=False)
        matingdata = newmatingmatxxxxxxxxxxx
        del newmatingmatxxxxxxxxxxx, axxxxxxxxxxxx, idfxxxxxxxxxxx
        matingdata.columns = ['SireID', 'DamID', 'SireIndex', 'DamIndex', 'PBTI', 'GEBV', 'MSV']
        matingdatafile = "{}/Mating data (Zygotic).csv".format(os.getcwd()) # output file
        matingdata.to_csv(matingdatafile, header = True, index=False) # prints Probdf of breeding top inds to file


    # Calculates RPTA
    elif (selstrat == 'RPTA'):
        RPTAdf = CalcRPTA(inputfile1, inputfile2, MSVMSC, index_wt, selinorthresh)
        matingdata = potentialparents(selstrat, RPTAdf, no_traits, selintmale, selintfemale, trait_names)
        GEBV = CalcGEBV(inputfile1, inputfile2, index_wt)
        if (no_traits == 1):
            matingdata.iloc[:, 5] = np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(0+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 6] = np.array(MSVMSC.iloc[:,0+2][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:,0+2][matingdata.iloc[:,3]])
            matingdata.iloc[:, 4] = matingdata.iloc[:, 5] + np.sqrt(matingdata.iloc[:, 6])*selinorthresh
        elif (no_traits > 1):
            matingdata.iloc[:, 5] = np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,2]])/2 + np.array(GEBV.iloc[:,(no_traits+2)][matingdata.iloc[:,3]])/2
            matingdata.iloc[:, 6] = np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,2]]) + np.array(MSVMSC.iloc[:, (no_traits*(no_traits-1))+2][matingdata.iloc[:,3]])
            matingdata.iloc[:, 4] = matingdata.iloc[:, 5] + np.sqrt(matingdata.iloc[:, 6])*selinorthresh

        idfxxxxxxxxxxx = np.unique(matingdata.iloc[:,3])
        newmatingmatxxxxxxxxxxx = pd.DataFrame(index=range(len(idfxxxxxxxxxxx)),columns=range(7))
        for mm in range(len(idfxxxxxxxxxxx)):
            axxxxxxxxxxxx = matingdata.loc[matingdata.iloc[:,3] == idfxxxxxxxxxxx[mm]]
            topsire = np.array(axxxxxxxxxxxx.iloc[:,2])
            newmatingmatxxxxxxxxxxx.iloc[mm, :] = axxxxxxxxxxxx.iloc[np.argmax(axxxxxxxxxxxx.iloc[:,4]),:]
            norepssire = Counter(newmatingmatxxxxxxxxxxx.iloc[:,2])
            for nr in range(len(topsire)):
                if (norepssire[topsire[nr]] <= 49):
                    newmatingmatxxxxxxxxxxx.iloc[mm, :] = np.array(axxxxxxxxxxxx[axxxxxxxxxxxx.iloc[:,2]==topsire[nr]])
                    break
        
        RPTAdffile = "{}/{}_gametic approach.csv".format(os.getcwd(), selstrat) # output file
        RPTAdf.to_csv(RPTAdffile, header = True, index=False)
        matingdata = newmatingmatxxxxxxxxxxx
        del newmatingmatxxxxxxxxxxx, axxxxxxxxxxxx, idfxxxxxxxxxxx
        matingdata.columns = ['SireID', 'DamID', 'SireIndex', 'DamIndex', 'RPTA', 'GEBV', 'MSV']
        matingdatafile = "{}/Mating data (Zygotic).csv".format(os.getcwd()) # output file
        matingdata.to_csv(matingdatafile, header = True, index=False) # prints RPTAdf of breeding top inds to file


    MSVMSCtemp = MSVMSC.iloc[:,2:MSVMSC.shape[1]]
    MSVMSCtemp = np.array(MSVMSCtemp.loc[matingdata.iloc[:,2], :]) + np.array(MSVMSCtemp.loc[matingdata.iloc[:,3], :])
    MSVMSCtemp = pd.DataFrame(MSVMSCtemp) #convert matrix to dataframe
    if (no_traits == 1):
        MSVMSCtemp.columns = trait_names
        MSVMSCtempfile = "{}/MSV_(zygotic).csv".format(os.getcwd()) # output file
    elif (no_traits > 1):
        tn = np.zeros((no_traits, no_traits), 'U20')
        tn = np.chararray(tn.shape, itemsize=20)
        for i in range(no_traits):
            for t in range(no_traits):
                if (i == t):
                    tn[i,t] = str(trait_names[i])
                elif (i != t):
                    tn[i,t] = "{}*{}".format(trait_names[i],trait_names[t])
        colnamesms = tn[np.tril_indices(no_traits)]
        colnamesms = np.concatenate((colnamesms, "Agg BV"), axis=None)
        MSVMSCtemp.columns = colnamesms
        MSVMSCtempfile = "{}/MSV_MSC_AggregateBV_(zygotic).csv".format(os.getcwd()) # output file
    MSVMSCtemp.insert(0, "Sire_ID", matingdata.iloc[:,0], True) # insert animal ID in the first col of the dataframe
    MSVMSCtemp.insert(1, "Dam_ID", matingdata.iloc[:,1], True) # insert animal ID in the second col of the dataframe
    MSVMSCtemp.to_csv(MSVMSCtempfile, header = True, index=False) # saves MSV and MSC and Aggre if mult traits
    return matingdata



# creates matrix for additive effects of aggregate genotype
def agggenmsezyg(inputfile1, inputfile2, Mlist, index_wt, ZygoDF):
    no_individuals1 = inputfile2.shape[0] #Number of individuals = nrow genotype matrix
    no_markers = inputfile1.shape[0] # no of markers = length of file
    MMfinal =  np.empty((no_individuals1, no_markers))
    for i in range(no_individuals1):
        tempmat1 = np.zeros((index_wt.size, no_markers))
        for tr in range(index_wt.size):
            tempmat1[tr,:] = Mlist[0][tr][i,:]
        MMfinal[i,:] = np.matmul(index_wt.transpose(), tempmat1) 
    MMfinal = pd.DataFrame(MMfinal)
    tempmatmale = np.array(MMfinal.loc[ZygoDF.iloc[:,2], :])
    tempmatfemale = np.array(MMfinal.loc[ZygoDF.iloc[:,3], :])
    return  tempmatmale, tempmatfemale 
        

# Derives similarity matrices using zygotic approach   
def SimMatZyg(inputfile1, inputfile2, index_wt, ZygoDF, CovMat, chrinterest, stdsim = True, progress = True):
    chrinterest = chrinterest
    if ('all' in chrinterest):
        chrinterest = 'all'
    elif ('none' in chrinterest):
        chrinterest = 'none'
    else:
        chrinterest = [int(i) for i in chrinterest]
        chrinterest = np.array(chrinterest)


    no_individuals = ZygoDF.shape[0] #Number of individuals = nrow genotype matrix
    no_chromosomes = max(inputfile1.iloc[:,1]) # no of chromosomes
    trait_names = inputfile1.columns[3:len(inputfile1.columns)] # traits names
    no_traits = trait_names.size
    
    psvlist = traitspecmatrices(inputfile1, inputfile2, no_traits) # list to store trait-specific matrices   
    
    # The following estimates similarity matrices
    #loop over each trait
    for i in range(no_traits):
        # create a matrix of zeros with rows and cols = no of inds
        mat = np.zeros((no_individuals, no_individuals)) #stores trait-specific covariance btw inds
        tempmat1xxxxxxxxx = pd.DataFrame(psvlist[0][i])
        for chr in range(no_chromosomes): #Loop over the number of chromosomes
            if progress:
                if (chr == 0):
                    print("Processing  ", trait_names[i])
            snpindex = np.array(inputfile1.iloc[:, 0][inputfile1.iloc[:, 1] == (1+chr)])-1 #marker index
            #matrix multiplication using BLAS
            tempmatmalexxxxxx = np.array(tempmat1xxxxxxxxx.loc[ZygoDF.iloc[:,2], snpindex])
            tempmatfemalexxxx = np.array(tempmat1xxxxxxxxx.loc[ZygoDF.iloc[:,3], snpindex])
            
            covtempxxxxxxxxxxxxxx = abs(dgemmMRM(tempmatmalexxxxxx, CovMat[0][chr])) + abs(dgemmMRM(tempmatfemalexxxx, CovMat[0][chr])) #((M*R)*M')
            mat = mat + covtempxxxxxxxxxxxxxx #sums up chr specific covariances
            del snpindex, tempmatmalexxxxxx, tempmatfemalexxxx
            if (type(chrinterest) == str):
                if (chrinterest == 'all'):
                    chrfile = "{}/Similarity matrix btw matepairs_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                    np.save(chrfile, covtempxxxxxxxxxxxxxx) # writes sim matrices for chromosomes of interest to file
            elif (chr+1 in chrinterest):
                chrfile = "{}/Similarity matrix btw matepairs_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                np.save(chrfile, covtempxxxxxxxxxxxxxx) # writes sim matrices for chromosomes of interest to file
            if stdsim:
                if (type(chrinterest) == str):
                    if (chrinterest == 'all'):
                        corrtempxxxxxxxxx = correlation_from_covariance(covtempxxxxxxxxxxxxxx)
                        chrfilec = "{}/Standardized similarity matrix btw matepairs_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                        np.save(chrfilec, corrtempxxxxxxxxx) # writes std sim matrices for chromosomes of interest to file
                        del corrtempxxxxxxxxx
                elif (chr+1 in chrinterest):
                    corrtempxxxxxxxxx = correlation_from_covariance(covtempxxxxxxxxxxxxxx)
                    chrfilec = "{}/Standardized similarity matrix btw matepairs_{}_chr_{}.npy".format(os.getcwd(), trait_names[i], chr+1)# output file
                    np.save(chrfilec, corrtempxxxxxxxxx) # writes std sim matrices for chromosomes of interest to file
                    del corrtempxxxxxxxxx
            del covtempxxxxxxxxxxxxxx
            #outputs progress in terms of trait name and chr processed
            if progress:
                progr(no_chromosomes, chr)
    del tempmat1xxxxxxxxx
    if (no_traits == 1):
        if stdsim:
            mat = correlation_from_covariance(mat)
    elif (no_traits > 1):
        if progress:
            print('Creating similarity matrix based on aggregate genotype')  
        tempmatmale, tempmatfemale = agggenmsezyg(inputfile1, inputfile2, psvlist, index_wt, ZygoDF)      
        mat = np.zeros((no_individuals, no_individuals)) # stores overall covariance btw inds
        # loop over chromososomes
        for chr in range(no_chromosomes):
            snpindex = np.array(inputfile1.iloc[:, 0][inputfile1.iloc[:, 1] == (1+chr)])-1
            covtempxxxxxxxxxxxxxx = abs(dgemmMRM(tempmatmale[:,snpindex],  CovMat[0][chr])) + abs(dgemmMRM(tempmatfemale[:,snpindex],  CovMat[0][chr]))
            mat = mat + covtempxxxxxxxxxxxxxx
            if progress: 
                progr(no_chromosomes, chr)

        if stdsim:
            mat = correlation_from_covariance(mat)
    return mat

# Optimizes mate allocation
def OptiMateAllocZyg(data, Similaritymat, maxalloc, progress = True):
    # Optimizing mate allocation
    no_individuals = data.shape[0]
    returns = data.iloc[:,4]
    muhat = np.linspace(start=min(returns), stop=max(returns), num=50)
    sln = np.zeros((len(muhat), len(returns)))
    opt_portfolios = np.zeros((len(muhat), 2))
    Similaritymat = matrix(Similaritymat)  # Covariance matrix
    q = matrix(np.zeros((len(returns), 1)))
    r_avg = np.zeros((1, len(returns)))
    r_avg[0,:] = returns   #selstrat values

    # Capturing inequality constraints Gx >= h
    # captures the constraints (r_avg'x >= muhat[i]) and (x >= lb and x <= ub)
    G = matrix(np.vstack((-r_avg, -1.0*np.eye(len(returns)) , np.eye(len(returns)))))
    L=np.array([-0.0]*len(returns)) #lower bound of 0.0
    # upper bound. to be provided from param file
    U=np.array([maxalloc]*len(returns))
    # Capturing equality constraint (Ax == b)
    A = matrix(1.0, (1,len(returns)))
    b = matrix(1.0)
    
    if progress:
        print('Optimizing portfolios')
        
    for i in range(len(muhat)):
        try:
            h=matrix(np.transpose(np.hstack((-muhat[i], L, U)).reshape((1,1+(2*len(returns))))))
            solvers.options['show_progress'] = False
            sol = solvers.qp(Similaritymat, q, G, h, A, b)
            sln[i,:] = np.array(sol['x']).transpose()
            opt_portfolios[i,0] = np.matmul(np.matmul(sln[i,:].transpose(), Similaritymat), sln[i,:])
            opt_portfolios[i,1] = np.dot(sln[i,:].transpose(), returns)
            
            if progress:
                progr(len(muhat), i)
        except Exception as err:
            exception_type = type(err).__name__
    if progress:
        if (maxalloc != 1):
            print("100", " ", flush=True)
        print('Mate allocation optimized')
    sln = sln[~np.all(sln == 0, axis=1)]
    slnfile = "{}/Optimum allocation of weights to individuals.csv".format(os.getcwd()) # output file
    sln = pd.DataFrame(sln)
    ID = data[['SireID', 'DamID']].agg('x'.join, axis=1)
    sln.columns = ID
    sln.to_csv(slnfile, header = True, index=False) # prints portfolio weights

    opt_portfolios = opt_portfolios[~np.all(opt_portfolios == 0, axis=1)]

    PlotEF(opt_portfolios, Similaritymat, returns, no_individuals)

    opt_portfolios = pd.DataFrame(opt_portfolios)
    opt_portfolios.columns = ['Risk', 'Returns']
    optfile = "{}/Expected risk and returns on portfolios.csv".format(os.getcwd()) # output file
    opt_portfolios.to_csv(optfile, header = True, index=False) # prints portfolio risks and returns
    return opt_portfolios, sln
