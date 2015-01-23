#!/usr/bin/env python 

#   file: dEcounters.py
#
#   Author:   Jia Liu    <liu.2053@osu.edu>
#   History:
#   Mar 03, 2014	 Include function getT0muFromiS() to calculate T0mu.
#   Feb 20, 2014	 Include rotation to T^0mu to convert 0 to local index
#   Feb 19, 2014     Disabled for now: Test mode: no initial flow: find 'jia test' for details.
#   Feb 18, 2014     Add dEdyViscousHydro() to calculate viscous hydro initial energy.
#   Feb 12, 2014     Add an argument iSFolder to dEdydphipSigmaFO(): let it find tables.
#   Feb 10, 2014     Add dEdyIdealHydro(): calculate initial ideal hydro part energy.
#   Feb 04, 2014     Add dEdyTotalInital(): calculate total initial energy.
#   Jan 31, 2014     dEdydphipSigmaFO() now use dE_dyptdptdphip table.
#   Jan 30, 2014     Add function dEdyTotalDecdat2() to calculate total energy
#				  flowing out of the freeze-out surface.
#   Jan 28, 2014     Fixed normalization bug of dEdydphipSigmaThermal.
#   Jan 24, 2014     Read the 'chosen_particles_backup.dat' to maintain
#    				 the correct matrix dimension. Commented dN calculation.

import sys,shutil
from os import path, stat, getcwd
import numpy as np

# some constants
Edec = 0.18  #Unit: GeV/fm^3
dxdy = 0.01
# table locations
rootDir = path.abspath('..')
tables_location = path.join(rootDir, 'tables')

def  dEdydphipSigmaThermal(dEdyd2rdphipFile, edFile, sfactor, \
		targetFolder, fileName):
	"""
	kick out all the elements outside the freeze-out surface and do
	integration on the transverse plane. Return total energy of the
	thermalization surface and dEdydphip_{\Sigma_{thermal}}.
	"""
	dEdyd2rdphip_data = np.loadtxt(dEdyd2rdphipFile)
	ed_data = np.loadtxt(edFile)
	phipTbl_file = path.join(tables_location, 'phip_gauss_table.dat')
	phipTbl = np.loadtxt(phipTbl_file) #100 points Gaussian points table for phip integration.
	phipGaussWeight = phipTbl[:, 1]

	# check dimensions
	ed_dataNum = ed_data.size
	if(ed_dataNum!=dEdyd2rdphip_data.shape[0]):
		print 'dEdydphipSigmaThermal error: Dimensions do not match!\n'
		print 'ed matrix elements: '+str(ed_dataNum), \
			', dEdyd2rdphip lines:' + str(dEdyd2rdphip_data.shape[0])
		sys.exit(-1)
	elif(phipGaussWeight.size!=dEdyd2rdphip_data.shape[1]):
		print 'dEdydphipSigmaThermal error: Dimensions do not match!\n'
		print 'phi_p Gaussian table length: '+str(phipGaussWeight.size),\
			', dEdyd2rdphip columns: ' + str(dEdyd2rdphip_data.shape[1])
		sys.exit(-1)

	# extract the row which are corresponding to cells outside the freeze-out 
	#surface
	ed_data = sfactor*ed_data
	ed_criteria = ed_data < Edec
	ed_criteria_rowNum = np.reshape(ed_criteria, (ed_dataNum)) #reshape is done row-wise
															# non-zero rows ed<Edec
	# integrate over the left cells
	dEdyd2rdphip_data = dEdyd2rdphip_data*sfactor
	dEdyd2rdphip_outside = dEdyd2rdphip_data[ed_criteria_rowNum, :]
	dEdydphip = dEdyd2rdphip_outside.sum(axis=0)*dxdy  # summing along each columnn
	dEdy = (dEdydphip*phipGaussWeight).sum()
	#print 'dEdyd2rdphip_outside', dEdyd2rdphip_outside.shape
	#print 'dEdydphip', dEdydphip.shape
	# save dEdydphip table
	savefileName = path.join(targetFolder, fileName)
	np.savetxt(savefileName, dEdydphip, fmt='%19.10e')
	return dEdy



def dEdydphipSigmaFO(iSFolder, iSDataFolder, targetFolder, fileName):
	"""
	calculate the total energy of the system from freeze-out surface (iS results)
	by summing over all thermal particles. Return the total transverse energy of
	freeze-out surface.
	"""
	# safty check, make sure iS has runned
	decdat2_file = path.join(iSDataFolder, 'decdat2.dat')
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0	
	# get table folders and files
	EOSTblFolder = path.join(iSFolder, 'EOS')
	GaussTblFolder = path.join(iSFolder, 'tables')

	chosenParticleTbl = np.loadtxt(path.join(EOSTblFolder,'chosen_particles_backup.dat'))#debug
	ptTbl = np.loadtxt(path.join(GaussTblFolder, 'pT_gauss_table.dat'))
	phipTbl = np.loadtxt(path.join(GaussTblFolder, 'phi_gauss_table.dat'))
	dEdyptdptdphipTbl = np.loadtxt(path.join(iSDataFolder, 'dE_ptdptdphidy.dat'))
	particleMassTbl = np.loadtxt(path.join(EOSTblFolder,'particles_mass.dat'))
	particleMassTbl = particleMassTbl[:,2]  #first column is particle index

	#check the dimension of data tables
	particleTotalNum = chosenParticleTbl.size
	phipTblLength = phipTbl.shape[0] #debug
	ptTblLength = ptTbl.shape[0]
	# delete gamma from dNptdptdphidy table
	dEdyptdptdphipTbl=dEdyptdptdphipTbl[phipTblLength::,:]

	# dimension check
	if(dEdyptdptdphipTbl.shape[0]!=particleTotalNum*phipTblLength):
		print 'Wrong phi_p table!\n '
		print 'Length of phip Gaussian table:'+str(phipTblLength), \
			 ', dEdyptdptdphip table length: '+str(dEdyptdptdphipTbl.shape[0]/particleTotalNum) 
		sys.exit(-1)
	elif(dEdyptdptdphipTbl.shape[1]!=ptTblLength):
		print 'Wrong pT table: \n '
		print 'length of pT Gaussian table: '+str(ptTblLength),\
			', dEdyptdptdphip table length: '+str(dEdyptdptdphipTbl.shape[1])
		sys.exit(-1)

	dEdy = 0
	dEdydphipTbl = np.zeros((phipTblLength), float) #for all particles
	#dNdy = 0
	#dNdydphipTbl = np.zeros((phipTblLength), float) #for all particles
	pT_mat = np.tile(ptTbl[:,0].transpose(), (phipTblLength, 1))
	pTGaussWeight_mat = np.tile(ptTbl[:,1].transpose(), (phipTblLength, 1)) # tile(M, (m,n)): repeat matrix for m,n times
	phipGaussWeight_array = phipTbl[:,1]

	for i in range(0, particleTotalNum):
		tempTbl = dEdyptdptdphipTbl[i*phipTblLength:phipTblLength*(i+1), :] #per hadron
		# do integration on pT: dNTbl.*pT.*mT.*weight(pT)
		particleMass = particleMassTbl[i]
		dEdydphipTemp = (pT_mat*tempTbl*pTGaussWeight_mat).sum(axis=1)  #summing along one row
		# dNdydphipTemp = (pT_mat*tempTbl*pTGaussWeight_mat).sum(axis=1)
		#print 'Size of dEdydphip table:'+str(dEdydphipTemp.shape)
		#print 'Size of dEdydphip table:'+str(dEdydphipTbl.shape)
		dEdydphipTbl = dEdydphipTbl+dEdydphipTemp
		# dNdydphipTbl = dNdydphipTbl+dNdydphipTemp

		# integrate over phip
		dEdy = dEdy+(dEdydphipTemp*phipGaussWeight_array).sum()
		# dNdy = dNdy+(dNdydphipTemp*phipGaussWeight_array).sum()
		
	# save to file
	savefileName=path.join(targetFolder, fileName)
	np.savetxt(savefileName, dEdydphipTbl, fmt='%19.10e')

	# savefileName=path.join(targetFolder, 'dNdydphip.dat')
	# np.savetxt(savefileName, dNdydphipTbl, fmt='%19.10e')

	return dEdy


def getTotaldEdy(dEdyd2rdphipFile, edFile, normFactor, iSFolder, iSdataFolder, \
        dEdphipOutName = 'dEdydphipThermal.dat',dEdphipHydroName = 'dEdydphipFO.dat'):
    """
    run two functions from dEcounters to get the total dEdy, also generate 
    dEdydphip tables for out surface and hydro surface. 
    """
    dEdy_therm = dEdydphipSigmaThermal(dEdyd2rdphipFile,edFile,\
        normFactor, iSdataFolder, dEdphipOutName)
    dEdy_fo = dEdydphipSigmaFO(iSFolder, iSdataFolder,\
        iSdataFolder, dEdphipHydroName)
    # total energy 
    totaldEdy = dEdy_therm + dEdy_fo
    return totaldEdy 



def dEdyTotalDecdat2(resultsFolder, fileName = 'decdat2.dat'):
	"""
	Calculate total energy from hydro freeze-out surface, i.e. file 'decdat2.dat'
	dEdy = \int T^{0\mu} d^3\sigma_\mu
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	# recombine T^0\mu
	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	ed_fo = decData_fo[:,6]
	pl_fo = decData_fo[:,11]
	ppi_fo = decData_fo[:,19]
	pi00_fo = decData_fo[:,13]
	pi01_fo = decData_fo[:,14]
	pi02_fo = decData_fo[:,15]
	pi11_fo = decData_fo[:, 16]
	pi12_fo = decData_fo[:, 17]
	pi22_fo = decData_fo[:, 18]

	T00_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*gamma_fo - pl_fo - ppi_fo \
		+ pi00_fo
	T01_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*ux_fo + pi01_fo		
	T02_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*uy_fo + pi02_fo 
	T0md3sigm = decData_fo[:,0]*(T00_fo*decData_fo[:,1] + T01_fo*decData_fo[:,2] \
		+ T02_fo*decData_fo[:,3])

	dEdy = T0md3sigm.sum()

	return dEdy


def dEdyTotalDecdat2_ideal(resultsFolder, fileName = 'decdat2.dat'):
	"""
	Calculate total energy from hydro freeze-out surface, i.e. file 'decdat2.dat'
	dEdy = \int T^{0\mu} d^3\sigma_\mu
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	# recombine T^0\mu
	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	ed_fo = decData_fo[:,6]
	pl_fo = decData_fo[:,11]

	T00_fo = (ed_fo+pl_fo)*gamma_fo*gamma_fo - pl_fo 
	T01_fo = (ed_fo+pl_fo)*gamma_fo*ux_fo 	
	T02_fo = (ed_fo+pl_fo)*gamma_fo*uy_fo  
	T0md3sigm = decData_fo[:,0]*(T00_fo*decData_fo[:,1] + T01_fo*decData_fo[:,2] \
		+ T02_fo*decData_fo[:,3])

	dEdy = T0md3sigm.sum()

	return dEdy

def dEdytotalConstSurface(resultsFolder, fileName = 'constTauSurface.dat',dxdy=0.01):
	"""
	Calculate total energy from hydro constant time surface, i.e. file 'decdat2.dat'
	dEdeta_s = \int T^{00} d^3\sigma_0
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	# recombine T^0\mu
	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	ed_fo = decData_fo[:,6]
	pl_fo = decData_fo[:,11]
	ppi_fo = decData_fo[:,19]
	pi00_fo = decData_fo[:,13]
	pi01_fo = decData_fo[:,14]
	pi02_fo = decData_fo[:,15]
	pi11_fo = decData_fo[:, 16]
	pi12_fo = decData_fo[:, 17]
	pi22_fo = decData_fo[:, 18]

	T00_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*gamma_fo - pl_fo - ppi_fo \
		+ pi00_fo
	T00d3sigm = decData_fo[:,0]*T00_fo*dxdy
 
	# include boost?
	global includeBoost
	if includeBoost == True:
		T01_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*ux_fo + pi01_fo		
		T02_fo = (ed_fo+pl_fo+ppi_fo)*gamma_fo*uy_fo + pi02_fo	
		T00_local = gamma_fo*T00_fo-ux_fo*T01_fo-uy_fo*T02_fo
		T00d3sigm = decData_fo[:,0]*T00_local*dxdy

	dEdy = T00d3sigm.sum()

	return dEdy

def radialFlowFO(resultsFolder, fileName = 'decdat2.dat'):
	"""
	Calculate radial flow on freeze-out surface
	formula: v_Sigma_fo=(\int_Sigma_fou^mu d^3sigma_mu v_perp)/(\int_Sigma_fo u^mu d^3sigma_mu)
	"""
	decdat2_file = path.join(resultsFolder, fileName)
	if(stat(decdat2_file).st_size==0):  # empty decdat2 file, hydro does not run
		return 0
	decData_fo = np.loadtxt(decdat2_file)

	vx_fo = decData_fo[:,4]
	vy_fo = decData_fo[:,5]
	gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
	ux_fo = decData_fo[:,4]*gamma_fo
	uy_fo = decData_fo[:,5]*gamma_fo
	vr_fo = np.sqrt(vx_fo**2+vy_fo**2)
	udsigma = decData_fo[:,0]*(gamma_fo*decData_fo[:,1]+\
		+ux_fo*decData_fo[:,2]+uy_fo*decData_fo[:,3])
	vr = sum(udsigma*vr_fo)/sum(udsigma)
	return vr




def dEdyTotalInital(dEdyphipFile, sfactor):
	"""
	calculate the total inital energy, by integrate the angular distribution dEdydphip
	"""
	try:
		dEdydphip_data = np.loadtxt(dEdyphipFile)
		phipTbl_file = path.join(tables_location, 'phip_gauss_table.dat')
		phipTbl = np.loadtxt(phipTbl_file) #100 points Gaussian points table for phip integration.
	except:
		print 'dEdyTotalInital: No such files!'
		sys.exit(-1)
	phipGaussWeight = phipTbl[:, 1]
	dEdydphip_data*=sfactor
	# dimension check
	if(dEdydphip_data.shape[0]!=phipGaussWeight.shape[0]):
		print 'Wrong phi_p table!\n '
		print 'Length of phip Gaussian table:'+str(phipGaussWeight.shape[0]), \
			 ', dEdydphip table length: '+str(dEdydphip_data.shape[0]) 
		sys.exit(-1)

	dEdy = (dEdydphip_data*phipGaussWeight).sum()
	return dEdy


def dEdyIdealHydro(InitialFolder, sfactor, tau, dxdy):
	"""
	calculate the intial total energy inside hydro surface 
	if viscous terms are dropped directly.
	"""
	try:
		ed_data  = np.loadtxt(path.join(InitialFolder, 'ed_profile_kln.dat'))
		pressure_data = np.loadtxt(path.join(InitialFolder, 'Pressure_kln.dat'))
		ux_data = np.loadtxt(path.join(InitialFolder, 'ux_profile_kln.dat'))
		uy_data = np.loadtxt(path.join(InitialFolder, 'uy_profile_kln.dat'))
	except:
		print 'dEdyIdealHydro: reading files error!'
		sys.exit(-1)
	ed_scaled = ed_data*sfactor
	ed_criteria = ed_scaled >= Edec
	gamma_init = np.sqrt(1 + ux_data**2 + uy_data**2)
#	gamma_init = np.ones(gamma_init.shape) # jia test
	T00_init = (ed_data + pressure_data)*gamma_init*gamma_init \
		 - pressure_data
	T00_init_outside = T00_init[ed_criteria]  #only count out-of-freezeout elements
	dEdy_init = (T00_init_outside.sum()).sum()*dxdy*sfactor*tau
	return dEdy_init

def dEdyViscousHydro(InitialFolder, sfactor, tau, dxdy):
	"""
	calculate the intial total energy inside hydro surface 
	including viscous terms.
	"""
	try:
		ed_data  = np.loadtxt(path.join(InitialFolder, 'ed_profile_kln.dat'))
		pressure_data = np.loadtxt(path.join(InitialFolder, 'Pressure_kln.dat'))
		ux_data = np.loadtxt(path.join(InitialFolder, 'ux_profile_kln.dat'))
		uy_data = np.loadtxt(path.join(InitialFolder, 'uy_profile_kln.dat'))
		ppi_data = np.loadtxt(path.join(InitialFolder, 'BulkPi_kln.dat'))
		pi00_data = np.loadtxt(path.join(InitialFolder, 'Pi00_kln.dat'))
	except:
		print 'dEdyViscousHydro: reading files error!'
		sys.exit(-1)
	ed_scaled = ed_data*sfactor
	ed_criteria = ed_scaled >= Edec
	gamma_init = np.sqrt(1 + ux_data**2 + uy_data**2)
	T00_init = (ed_data + pressure_data + ppi_data)*gamma_init*gamma_init \
		 - pressure_data - ppi_data + pi00_data
	T00_init_outside = T00_init[ed_criteria]  #only count out-of-freezeout elements
	dEdy_init = (T00_init_outside.sum()).sum()*dxdy*sfactor*tau
	return dEdy_init

def getT0muFromiS(iSEOSFolder, iSDataFolder):
	"""
	Find T0mu of Cooper-Frye by summing over all thermal particles.
	"""
	try:
		chosenParticleTbl = np.loadtxt(path.join(iSEOSFolder, 'chosen_particles.dat'))
	except:
		print 'getT0muFromiS: cannot load chosen_particle table!'
	T0mu = np.zeros((4,1))
	i=0
	for particle_idx in chosenParticleTbl:
		thermalParticleETFile = 'thermal_%d_ET_integrated_vndata.dat' %particle_idx
		try:
			thermalParticleTbl = np.loadtxt(path.join(iSDataFolder, thermalParticleETFile))
		except:
			print 'getT0muFromiS: cannot load ET file for: '+str(particle_idx)
		for j in range(0,4):
			T0mu[j] += thermalParticleTbl[0, j+2]
		i+=1
	return T0mu[:]


def doubleProjectPimn(decdat2_line):
    """
    Applying double projection operator on hydro dumped pi^munu.
    Return the line with replaced pi^munu
    """
    gmn = np.mat([[1., 0., 0., 0.], \
        [0., -1, 0., 0.],   \
        [0., 0., -1., 0.],  \
        [0., 0., 0., -1.]]) # metric
    # get umu rank-1 matrix
    vx_fo = decdat2_line[4]
    vy_fo = decdat2_line[5]
    gamma_fo = 1.0/np.sqrt(1.0 - vx_fo**2- vy_fo**2 )
    ux_fo = decdat2_line[4]*gamma_fo
    uy_fo = decdat2_line[5]*gamma_fo    
    uz_fo = 0
    umu = np.mat([gamma_fo, ux_fo, uy_fo, uz_fo])
    # get pi^munu matrix
    pimn = np.zeros([1, 7])
    for i in range(1,7):
        pimn[0,i-1] = decdat2_line[i+12] # pi00, pi01, pi02, pi11, pi12, pi22
    pimn[0,-1] = decdat2_line[12] # pi33
    pimn_input = np.mat([[pimn[0, 0], pimn[0, 1], pimn[0, 2], 0.],    \
    [pimn[0, 1], pimn[0, 3], pimn[0, 4], 0.], \
    [pimn[0, 2], pimn[0, 4], pimn[0, 5], 0], \
    [0., 0., 0., pimn[0, 6]] ])

    # apply double projection operator on pi^munu
    pimn_ap = (gmn*gmn - umu.T*(umu*gmn))*(((gmn*gmn - umu.T*(umu*gmn))*pimn_input).T) \
        -1.0/3.0*(gmn - umu.T*umu)*((np.trace(pimn_input*gmn)-(umu*gmn)*pimn_input*(gmn*umu.T)))[0,0]
    # update pi^munu the decdat2_oneLine
    decdat2_line[12] = pimn_ap[3,3]
    decdat2_line[13] = pimn_ap[0,0]
    decdat2_line[14] = pimn_ap[0,1]
    decdat2_line[15] = pimn_ap[0,2]
    decdat2_line[16] = pimn_ap[1,1]
    decdat2_line[17] = pimn_ap[1,2]
    decdat2_line[18] = pimn_ap[2,2]
    return decdat2_line


def preProcessDecdat2File(source_folder):
    """
    Pre-process input decdat2 file, replacing the
    original decdat2.dat. Backup the original one.
    """
    decdat2_fileName = 'decdat2.dat' 
    decdat2_file = path.join(source_folder, decdat2_fileName)
    try:
        decdat2_input = np.loadtxt(decdat2_file)
    except:
        print 'preProcessPimn: cannot read file: '+decdat2_file
        sys.exit(-1)

    # process file
    if decdat2_input.ndim == 1: # only one line in decdat2.dat
        decdat2_oneLine = decdat2_input
        decdat2_output = doubleProjectPimn(decdat2_oneLine)
    else:
        decdat2_output = np.zeros((1, 20)) # prepare a zero line to accept the following data
        decdat2_file_length = decdat2_input.shape[0]
        for i in range(0, decdat2_file_length):
            decdat2_oneLine = decdat2_input[i,:]
            decdat2_oneLine_new = doubleProjectPimn(decdat2_oneLine)
            decdat2_output = np.concatenate((decdat2_output,np.mat(decdat2_oneLine_new)))
        decdat2_output = np.delete(decdat2_output, 0, axis=0) # delete the zero data in line 1

    # back up the old decdat2.dat and save the new one
    decdat2_backupFileName = 'decdat2_backup.dat'
    decdat2_backupFile = path.join(source_folder, decdat2_backupFileName)
    shutil.move(decdat2_file, decdat2_backupFile)
    if decdat2_input.ndim == 1:
        np.savetxt(decdat2_file, decdat2_output[None], fmt='%14.6E')
    else:
        np.savetxt(decdat2_file, decdat2_output, fmt='%14.6E')


def calculateWn(folder, dEdphipOutName = 'dEdydphipThermal.dat',
    dEdphipHydroName = 'dEdydphipFO.dat', outputFileName = 'wn_integrated_vndata.dat'):
    """
    This function calculate energy flow anisotropy omega_n (Wn) from the energy azimuthal distribution
    at the hydro surface and the out surface.
    It is the python version of part of the matlab script "" which calculates the same quantity.
    The output Format is has 10 rows and 5 columns for order 1 to order 9
    """
    try:
        # load phip table
        phip_th = np.loadtxt(path.join(tables_location, 'phip_gauss_table.dat'))
        phip_fo = np.loadtxt(path.join(tables_location, 'phi_gauss_table.dat'))
        # load energy distribution
        dEdydphip_fo_data = np.loadtxt(path.join(folder, dEdphipHydroName))
        dEdydphip_th_data = np.loadtxt(path.join(folder, dEdphipOutName))
    except:
        print "calculateWn: no necessary files. Check iS/results Folder!\n"
        return False

    try:
        # calculate total energy
        wn_th_denominator = np.sum(dEdydphip_th_data*phip_th[::,1])
        wn_fo_denominator = np.sum(dEdydphip_fo_data*phip_fo[::,1])
    except:
        print "calculateWn: dimensions do not match. Check phi and phip tables!\n"
        return False

    # pre-allocate space
    wn_data = np.zeros((10, 6))

    # loop over all orders
    for iorder in range(0,10):
        # numerators for energy at out surface
        wn_th_numerator_real = np.sum(dEdydphip_th_data*np.cos(iorder*phip_th[:,0])*phip_th[:,1])
        wn_th_numerator_img  = np.sum(dEdydphip_th_data*np.sin(iorder*phip_th[:,0])*phip_th[:,1])
        # numerators for energy at hydro surface
        wn_fo_numerator_real = np.sum(dEdydphip_fo_data*np.cos(iorder*phip_fo[:,0])*phip_fo[:,1])
        wn_fo_numerator_img  = np.sum(dEdydphip_fo_data*np.sin(iorder*phip_fo[:,0])*phip_fo[:,1])
        # calculate omega n
        wn_real = (wn_th_numerator_real+wn_fo_numerator_real)/(wn_th_denominator+wn_fo_denominator+1e-18)
        wn_img  = (wn_th_numerator_img + wn_fo_numerator_img)/(wn_th_denominator+wn_fo_denominator+1e-18)
        wn_data[iorder, :] = np.array([iorder, (wn_th_numerator_real+wn_fo_numerator_real),
                                               (wn_th_numerator_img +wn_fo_numerator_img),
                                               wn_real, wn_img, np.sqrt(wn_real**2.0+wn_img**2.0)])
    # save file to folder
    savefileName = path.join(folder, outputFileName)
    np.savetxt(savefileName, wn_data,fmt='%19.8e')
    return True



