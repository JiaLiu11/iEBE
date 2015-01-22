#!/usr/bin/env python2
"""
    This module finds or updates the scaling factor for hydro run from sfactorDB.db.

"""
import sqlite3, sys

def sfactorFinder(db_name, values_list):
	"""
	This function finds the sfactor correspoinding to specific 
	switching time tau_s, shear viscosity eta/s (eta_s), and decoupling 
	energy density e_dec
	"""
	# create connection to database
	dbCon = sqlite3.connect(db_name)
	conCursor = dbCon.cursor()

	# refine information
	names_list = ['tau_s', 'eta_s', 'e_dec', 'model', 'has_fs']
	values_list_str = [str(i) for i in values_list]
	nameValue_list = zip(names_list[0:3], values_list_str[0:3])

	# generate query
	select_clause = "select sfactor from sfactorTable "
	where_clause = " where "+ " and ".join(map("=".join, nameValue_list[0:3]))+ \
					" and model = '%s' and has_fs = '%s' "%(values_list_str[3], 
															values_list_str[4])
	sql_query = select_clause + where_clause
	sfactor = conCursor.execute(sql_query).fetchall()

	# close connection
	dbCon.close()
	
	# check result
	if len(sfactor)==0:
		print "sfactorUtility: No sfactor found, please check the input!"
		print "parameter names: ", names_list
		print "values:          ", values_list_str, "\n"
		sys.exit(-1)		
	elif len(sfactor)>1: 
		print "sfactorUtility: More than one sfactor, please check the input!"
		print "parameter names: ", names_list
		print "values:          ", values_list_str,"\n"
		sys.exit(-1)

	return sfactor[0][0]

