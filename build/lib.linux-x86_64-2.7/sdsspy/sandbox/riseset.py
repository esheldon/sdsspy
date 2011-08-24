################################################################################
#
# RiseSet.py - This module provides functions to easily get the rise and
#				set time of a given RA and Dec, given a geographic latitude
#				and longitude, and date.
#
# Copyright (C) 2011 Adrian Price-Whelan
# 
################################################################################

"""
todo: 
"""

__author__ = 'Adrian Price-Whelan <adrn@nyu.edu>'

# Standard library dependencies (e.g. sys, os)
import datetime
from math import cos, sin, acos, asin, degrees, radians, tan, modf

from apwlib.Geometry import Angle
from apwlib.astro.convert import *

def rise(ra, dec, lat, lon, date=datetime.datetime.now(), raUnits='DEGREES', decUnits='DEGREES', latUnits='DEGREES', lonUnits='DEGREES', lonDirection='W'):
	""" 
	This function calculates the rise time of an object at a specified RA, 
	Dec, and Geographic Latitude. This is adapted from Practical Astronomy
	with your Calculator, 3rd Edition, page 53.
	
	Parameters
	----------
	ra : float, int, apwlib.Angle
		A Right Ascension - if float or int, pay attention to raUnits.
	dec : float, int, apwlib.Angle
		A Declination - if float or int, pay attention to decUnits.
	lat : float, int, apwlib.Angle
		A Geographic Latitude - if float or int, pay attention to latUnits.
	lon : float, int, apwlib.Angle
		A Geographic Longitude - if float or int, pay attention to lonUnits and lonDirection.
	date : datetime.datetime or string
		Can either be a datetime object or a string of an ISO 8601 date (e.g. 2010-01-24), 
		from which the date will be taken, and used to compute the rise time of the input RA 
		and Dec on that date. 
	raUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	decUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	latUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	lonUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	lonDirection : string
		Default is longitude WEST, 'W', but you can specify EAST by passing 'E'.
		
	Returns
	-------
	rise : datetime.datetime
		Returns a datetime.datetime object corresponding to the rise time 
		in UT of the input RA and Dec at the specified Latitude.
	 
	Examples
	--------
	>>> alpha = Angle()
	>>> alpha.hms = "23:39:20"
	>>> delta = Angle()
	>>> delta.deg = sex2dec(21, 42, 0)
	>>> date = datetime.datetime(1980, 8, 24)
	>>> lat = Angle()
	>>> lat.deg = 30
	>>> lon = Angle()
	>>> lon.deg = 64
	>>> rise(alpha, delta, lat, lon, date, lonDirection='E')
	1980-08-24 14:18:08.805547
	
	>>> rise(354.8333333, 21.7, 30, 64, '1980-08-24', lonDirection='E')
	1980-08-24 14:18:08.805539
	
	"""
	if not isinstance(ra, Angle):
		alpha = Angle()
		if raUnits.upper() == 'HOURS':
			alpha.hours = ra			
		elif raUnits.upper() == 'DEGREES':
			alpha.deg = ra
		else:
			raise AssertionError('raUnits must be either HOURS or DEGREES')
	else:
		alpha = ra
	
	if not isinstance(dec, Angle):
		delta = Angle()
		if decUnits.upper() == 'HOURS':
			delta.hours = dec
		elif decUnits.upper() == 'DEGREES':
			delta.deg = dec
		else:
			raise AssertionError('decUnits must be either HOURS or DEGREES')
	else:
		delta = dec
	
	if not isinstance(lat, Angle):
		phi = Angle()
		if latUnits.upper() == 'HOURS':
			phi.hours = lat
		elif latUnits.upper() == 'DEGREES':
			phi.deg = lat
		else:
			raise AssertionError('latUnits must be either HOURS or DEGREES')
	else:
		phi = lat
	
	if not isinstance(lon, Angle):
		lambdaa = Angle()
		if lonUnits.upper() == 'HOURS':
			lambdaa.hours = lon
		elif lonUnits.upper() == 'DEGREES':
			lambdaa.deg = lon
		else:
			raise AssertionError('lonUnits must be either HOURS or DEGREES')
		
		if lonDirection == 'W':
			lambdaa.deg = -1.0 * lambdaa.deg
	else:
		lambdaa = lon
	
	if not (isinstance(date, datetime.datetime) or isinstance(date, datetime.date)):
		matchStr = re.compile("(\d\d\d\d)(\-|\s*)(\d\d)(\-|\s*)(\d\d)")
		dateSearch = matchStr.search(date)
		try: 
			y, sp1, m, sp2, d = dateSearch.groups()
		except AttributeError:
			raise ValueError("Date not formatted properly.")
		datetimeObj = datetime.datetime(int(y), int(m), int(d))
		
	else:
		datetimeObj = date
	
	# Compute the azimuths for rising and setting
	Ar = Angle()
	Ar.rad = acos(sin(delta.rad)/cos(phi.rad))
	
	HA = Angle()
	HA.deg = (degrees(acos(-tan(phi.rad)*tan(delta.rad)))) % 360.0
	
	LSTr = (24 + alpha.hours - HA.hours) % 24.0
	GMST = LSTr - lambdaa.hours
	
	h, m, sFloat = dec2sex(GMST)
	s2, s = modf(sFloat)
	micros = s2*10**6
	
	gmstDatetime = datetime.datetime(datetimeObj.year,datetimeObj.month,datetimeObj.day, h, m, int(s), int(micros))
	
	return gmst2utcDatetime(gmstDatetime)

def set(ra, dec, lat, lon, date=datetime.datetime.now(), raUnits='DEGREES', decUnits='DEGREES', latUnits='DEGREES', lonUnits='DEGREES', lonDirection='W'):
	""" 
	This function calculates the set time of an object at a specified RA, 
	Dec, and Geographic Latitude. This is adapted from Practical Astronomy
	with your Calculator, 3rd Edition, page 53.
	
	Parameters
	----------
	ra : float, int, apwlib.Angle
		A Right Ascension - if float or int, pay attention to raUnits.
	dec : float, int, apwlib.Angle
		A Declination - if float or int, pay attention to decUnits.
	lat : float, int, apwlib.Angle
		A Geographic Latitude - if float or int, pay attention to latUnits.
	lon : float, int, apwlib.Angle
		A Geographic Longitude - if float or int, pay attention to lonUnits and lonDirection.
	date : datetime.datetime or string
		Can either be a datetime object or a string of an ISO 8601 date (e.g. 2010-01-24), 
		from which the date will be taken, and used to compute the rise time of the input RA 
		and Dec on that date. 
	raUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	decUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	latUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	lonUnits : string
		Accepted values are DEGREES, HOURS. Default is DEGREES.
	lonDirection : string
		Default is longitude WEST, 'W', but you can specify EAST by passing 'E'.
		
	Returns
	-------
	rise : datetime.datetime
		Returns a datetime.datetime object corresponding to the set time 
		in UT of the input RA and Dec at the specified Latitude.
	 
	Examples
	--------
	>>> alpha = Angle()
	>>> alpha.hms = "23:39:20"
	>>> delta = Angle()
	>>> delta.deg = sex2dec(21, 42, 0)
	>>> date = datetime.datetime(1980, 8, 24)
	>>> lat = Angle()
	>>> lat.deg = 30
	>>> lon = Angle()
	>>> lon.deg = 64
	>>> set(alpha, delta, lat, lon, date, lonDirection='E')
	1980-08-24 04:06:05.039232
	
	>>> set(354.8333333, 21.7, 30, 64, '1980-08-24', lonDirection='E')
	1980-08-24 04:06:05.039224
	
	"""
	if not isinstance(ra, Angle):
		alpha = Angle()
		if raUnits.upper() == 'HOURS':
			alpha.hours = ra			
		elif raUnits.upper() == 'DEGREES':
			alpha.deg = ra
		else:
			raise AssertionError('raUnits must be either HOURS or DEGREES')
	else:
		alpha = ra
	
	if not isinstance(dec, Angle):
		delta = Angle()
		if decUnits.upper() == 'HOURS':
			delta.hours = dec
		elif decUnits.upper() == 'DEGREES':
			delta.deg = dec
		else:
			raise AssertionError('decUnits must be either HOURS or DEGREES')
	else:
		delta = dec
	
	if not isinstance(lat, Angle):
		phi = Angle()
		if latUnits.upper() == 'HOURS':
			phi.hours = lat
		elif latUnits.upper() == 'DEGREES':
			phi.deg = lat
		else:
			raise AssertionError('latUnits must be either HOURS or DEGREES')
	else:
		phi = lat
	
	if not isinstance(lon, Angle):
		lambdaa = Angle()
		if lonUnits.upper() == 'HOURS':
			lambdaa.hours = lon
		elif lonUnits.upper() == 'DEGREES':
			lambdaa.deg = lon
		else:
			raise AssertionError('lonUnits must be either HOURS or DEGREES')
		
		if lonDirection == 'W':
			lambdaa.deg = -1.0 * lambdaa.deg
	else:
		lambdaa = lon
	
	if not (isinstance(date, datetime.datetime) or isinstance(date, datetime.date)):
		matchStr = re.compile("(\d\d\d\d)(\-|\s*)(\d\d)(\-|\s*)(\d\d)")
		dateSearch = matchStr.search(date)
		try: 
			y, sp1, m, sp2, d = dateSearch.groups()
		except AttributeError:
			raise ValueError("Date not formatted properly.")
		datetimeObj = datetime.datetime(int(y), int(m), int(d))
		
	else:
		datetimeObj = date
	
	# Compute the azimuths for rising and setting
	Ar = Angle()
	Ar.rad = acos(sin(delta.rad)/cos(phi.rad))
	
	HA = Angle()
	HA.deg = (degrees(acos(-tan(phi.rad)*tan(delta.rad)))) % 360.0
	
	LSTs = (alpha.hours + HA.hours) % 24.0
	GMST = LSTs - lambdaa.hours
	
	h, m, sFloat = dec2sex(GMST)
	s2, s = modf(sFloat)
	micros = s2*10**6
	
	gmstDatetime = datetime.datetime(datetimeObj.year,datetimeObj.month,datetimeObj.day, h, m, int(s), int(micros))
	
	return gmst2utcDatetime(gmstDatetime)

if __name__ == '__main__':
	pass
	