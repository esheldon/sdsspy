################################################################################
#
# convert.py - A Python module that contains unit conversions and 
#				coordinate conversions for astronomical data
#				processing and analysis
#
# Copyright (C) 2011 Adrian Price-Whelan
#
################################################################################

__author__ = 'Adrian Price-Whelan <adrn@nyu.edu>'

# Standard library dependencies
import datetime
import math
import re

import apwlib.astro.geometry as geometry

#
# Assumed Constants:

# T is the number of Julian days from 2000 January 1.5, but we use that epoch for RA/Dec nowadays..
epochOffsetT = 0.0
ee = 23.0 + 26.0/60.0 + 21.45/3600.0 - 46.815/3600.0*epochOffsetT - 0.0006/3600.0*epochOffsetT**2 + 0.00181/3600.0*epochOffsetT**3

#
# Formatting conversions:

def dec2sex(deci):
	""" 
	Converts a Decimal number (in hours or degrees) to Sexagesimal.
	
	Parameters
	----------
	deci : float
		A decimal number to be converted to Sexagismal.
	
	Returns
	-------
	hd : int
		hours or degrees
	m : int
		minutes or arcminutes
	s : float
		seconds or arcseconds
	
	"""
	(hfrac, hd) = math.modf(deci)
	(min_frac, m) = math.modf(hfrac * 60)
	s = min_frac * 60.
	return (int(hd), int(m), s)

def sex2dec(hd, min, sec):
	""" 
	Converts a Sexagesimal number to a Decimal number.
	
	Parameters
	----------
	hd : int
		hours or degrees
	m : int
		minutes or arcminutes
	s : float
		seconds or arcseconds
	
	Returns
	-------
	hd : float
		A decimal number
	
	"""
	return float(hd) + min/60.0 + sec/3600.0

def datetime2decHours(time):
	""" 
	Converts a datetime.time or datetime.datetime object into decimal time.
	
	Parameters
	----------
	time : datetime.time or datetime.datetime
	
	Returns
	-------
	decTime : float
		A decimal number representing the input time
	
	"""
	return time.hour + time.minute/60.0 + time.second/3600.0 + time.microsecond/3600000000.0

def string2hours(string):
	""" 
	Converts a Sexagesimal string to a decimal hour
	
	Parameters
	----------
	string : str
		A string representing a sexagismal time, in the form HH:MM:SS.SS.
		
	Returns
	-------
	decHours : float
		The converted sexagismal time to decimal time.
	
	"""
	string = string.strip() # trim leading/trailing whitespace

	div = '[|:|\s]' # allowable field delimiters "|", ":", whitespace
	
	# (optional +/-) degrees + min + decimal seconds
	sdec = '([+-]?)(\d{1,3})' + div + '(\d{1,2})' + div + '(\d{1,2}\.?\d+?)'
	
	co_re= re.compile(sdec)
	co_search= co_re.search(string)
	
	if co_search is None:
		raise ValueError("Invalid input string: '%s'" % string)
	elems = co_search.groups()

	plusminus = elems[0]
	hours = float(elems[1])
	minutes = float(elems[2])
	seconds = float(elems[3])
	
	# Check for nonsense values
	if hours > 24.0:
		raise ValueError("Hour value must be < 24.")
	if minutes >= 60.0:
		raise ValueError("Minute value must be < 60.")
	if seconds >= 60.0:
		raise ValueError("Second value must be < 60.")
		
	# Convert dec
	decHours = hours + minutes/60.0 + seconds/3600.0
	
	if plusminus is "-":
		decHours = -1.0 * decHours
	
	return decHours
	
def string2deg(string):
	""" 
	Converts a Sexagesimal string to a decimal degrees
	
	Parameters
	----------
	string : str
		A string representing a sexagismal time, in the form [+|-]DD:MM:SS.SS
		
	Returns
	-------
	decHours : float
		The converted sexagismal degrees to decimal degrees
	
	"""
	string = string.strip() # trim leading/trailing whitespace

	div = '[|:|\s]' # allowable field delimiters "|", ":", whitespace
	
	# (optional +/-) degrees + min + decimal seconds
	sdec = '([+-]?)(\d{1,3})' + div + '(\d{1,2})' + div + '(\d{1,2}\.?\d+?)'
	
	co_re= re.compile(sdec)
	co_search= co_re.search(string)
	
	if co_search is None:
		raise ValueError("Invalid input string: %s" % string)
	elems = co_search.groups()

	plusminus = elems[0]
	degrees = float(elems[1])
	arcminutes = float(elems[2])
	arcseconds = float(elems[3])
	
	# Check for nonsense values
	if degrees > 90.0:
		raise ValueError("Degree value must be <= 90.")
	if arcminutes >= 60.0:
		raise ValueError("Arcminute value must be < 60.")
	if arcseconds >= 60.0:
		raise ValueError("Arcsecond value must be < 60 (was %f)." % arcseconds)
		
	# Convert dec
	decDegrees = degrees + arcminutes/60.0 + arcseconds/3600.0
	
	if plusminus is "-":
		decDegrees = -1.0 * decDegrees
	
	return decDegrees

#
# Date / Time Conversions:

def ymd2jd(year, month, day):
	""" 
	Converts a year, month, and day to a Julian Date.
	This function uses an algorithm from the book "Practical Astronomy with your 
	Calculator" by Peter Duffet-Smith (Page 7)
	
	Parameters
	----------
	year : int
		A Gregorian year
	month : int
		A Gregorian month
	day : int
		A Gregorian day
	
	Returns
	-------
	jd : float
		A Julian Date computed from the input year, month, and day.
	
	"""
	if month == 1 or month == 2:
		yprime = year - 1
		mprime = month + 12
	else:
		yprime = year
		mprime = month
	
	if year > 1582 or (year == 1582 and (month >= 10 and day >= 15)):
		A = int(yprime / 100)
		B = 2 - A + int(A/4.0)
	else:
		B = 0	
		
	if yprime < 0:
		C = int((365.25 * yprime) - 0.75)
	else:
		C = int(365.25 * yprime)
	
	D = int(30.6001 * (mprime + 1))
	
	return B + C + D + day + 1720994.5

def utcDatetime2gmst(datetimeObj):
	""" 
	Converts a Python datetime object with UTC time to Greenwich Mean Sidereal Time.
	This function uses an algorithm from the book "Practical Astronomy with your 
	Calculator" by Peter Duffet-Smith (Page 17)
	
	Parameters
	----------
	datetimeObj : datetime.datetime
		A Python datetime.datetime object
	
	Returns
	-------
	< > : datetime.datetime
		A Python datetime.datetime object corresponding to the Greenwich Mean 
		Sidereal Time of the input datetime.datetime object. 
	
	"""
	jd = ymd2jd(datetimeObj.year, datetimeObj.month, datetimeObj.day)
	
	S = jd - 2451545.0
	T = S / 36525.0
	T0 = 6.697374558 + (2400.051336 * T) + (0.000025862 * T**2)
	T0 = T0 % 24
	
	UT = datetime2decHours(datetimeObj.time()) * 1.002737909
	T0 += UT
	
	GST = T0 % 24

	h,m,s = dec2sex(GST)
	return datetime.datetime(year=datetimeObj.year, month=datetimeObj.month, day=datetimeObj.day, hour=h, minute=m, second=int(s), microsecond=int((s-int(s))*10**6))
	
def gmst2utcDatetime(datetimeObj):
	""" 
	Converts a Python datetime object representing a Greenwich Mean Sidereal Time 
	to UTC time. This function uses an algorithm from the book "Practical Astronomy
	with your Calculator" by Peter Duffet-Smith (Page 18)
	
	Parameters
	----------
	datetimeObj : datetime.datetime
		A Python datetime.datetime object
	
	Returns
	-------
	< > : datetime.datetime
		A Python datetime.datetime object corresponding to UTC time of the input 
		Greenwich Mean Sidereal Time datetime.datetime object. 
	
	"""
	jd = ymd2jd(datetimeObj.year, datetimeObj.month, datetimeObj.day)
	
	S = jd - 2451545.0
	T = S / 36525.0
	T0 = 6.697374558 + (2400.051336*T) + (0.000025862*T**2)
	T0 = T0 % 24
	
	GST = (datetime2decHours(datetimeObj.time()) - T0) % 24
	UT = GST * 0.9972695663
	
	h,m,s = dec2sex(UT)
	return datetime.datetime(year=datetimeObj.year, month=datetimeObj.month, day=datetimeObj.day, hour=h, minute=m, second=int(s), microsecond=int((s-int(s))*10**6))

def mjd2jd(mjd):
	""" 
	Converts a Modified Julian Date to Julian Date
	
	Parameters
	----------
	mjd : float (any numeric type)
		A Modified Julian Date
	
	Returns
	-------
	jd : float
		The Julian Date calculated from the input Modified Julian Date
	 
	Examples
	--------
	>>> mjd2jd(55580.90429)
	2455581.40429
	
	"""
	return mjd + 2400000.5

def jd2mjd(jd):
	""" 
	Converts a Julian Date to a Modified Julian Date
	
	Parameters
	----------
	jd : float (any numeric type)
		A julian date
	
	Returns
	-------
	mjd : float
		The Modified Julian Date (MJD) calculated from the input Julian Date
	 
	Examples
	--------
	>>> jd2mjd(2455581.40429)
	55580.90429
	
	"""
	return float(jd - 2400000.5)

def mjd2sdssjd(mjd):
	""" 
	Converts a Modified Julian Date to a SDSS Julian Date, which is
		used at Apache Point Observatory (APO).
	
	Parameters
	----------
	mjd : float (any numeric type)
		A Modified Julian Date
	
	Returns
	-------
	mjd : float
		The SDSS Julian Date calculated from the input Modified Julian Date
	
	Notes
	-----
	- The SDSS Julian Date is a convenience used at APO to prevent MJD from rolling
		over during a night's observation. 
	
	Examples
	--------
	>>> mjd2sdssjd(55580.90429)
	55581.20429
	
	"""
	return mjd + 0.3
	
def jd2sdssjd(jd):
	""" 
	Converts a Julian Date to a SDSS Julian Date, which is
		used at Apache Point Observatory (APO).
	
	Parameters
	----------
	jd : float (any numeric type)
		A Julian Date
	
	Returns
	-------
	mjd : float
		The SDSS Julian Date calculated from the input Modified Julian Date
	
	Notes
	-----
	- The SDSS Julian Date is a convenience used at APO to prevent MJD from rolling
		over during a night's observation. 
	
	Examples
	--------
	>>> jd2sdssjd(2455581.40429)
	55581.20429
	
	"""
	mjd = jd2mjd(jd)
	return mjd2sdssjd(mjd)

def sdssjd2mjd(sdssjd):
	""" 
	Converts a SDSS Julian Date to a Modified Julian Date 
	
	Parameters
	----------
	sdssjd : float (any numeric type)
		The SDSS Julian Date
	
	Returns
	-------
	mjd : float
		The Modified Julian Date calculated from the input SDSS Julian Date
	
	Notes
	-----
	- The SDSS Julian Date is a convenience used at APO to prevent MJD from rolling
		over during a night's observation. 
	
	Examples
	--------
	>>> sdssjd2mjd(55581.20429)
	55580.90429
	
	"""
	return sdssjd - 0.3
	
def sdssjd2jd(sdssjd):
	""" 
	Converts a SDSS Julian Date to a Julian Date 
	
	Parameters
	----------
	sdssjd : float (any numeric type)
		The SDSS Julian Date
	
	Returns
	-------
	jd : float
		The Julian Date calculated from the input SDSS Julian Date
	
	Notes
	-----
	- The SDSS Julian Date is a convenience used at APO to prevent MJD from rolling
		over during a night's observation. 
	
	Examples
	--------
	>>> sdssjd2jd(55581.20429)
	2455581.40429
	
	"""
	mjd = sdssjd - 0.3
	return mjd2jd(mjd)

def jd2datetime(jd):
	""" 
	Converts a Julian Date to a Python datetime object. The resulting time is in UTC.
	
	Parameters
	----------
	jd : float (any numeric type)
		A Julian Date
	
	Returns
	-------
	< > : datetime.datetime
		A Python datetime.datetime object calculated using an algorithm from the book
		"Practical Astronomy with your Calculator" by Peter Duffet-Smith (Page 8)
	
	Examples
	--------
	>>> jd2datetime(2455581.40429)
	2011-01-19 21:42:10.000010
	
	"""
	
	jd = jd + 0.5
	(jd_frac, jd_int) = math.modf(jd)
	
	if jd_int > 2299160:
		A = int((jd_int - 1867216.25)/36524.25)
		B = jd_int + 1 + A - int(A/4)
	else:
		B = jd_int
		
	C = B + 1524
	D = int((C - 122.1) / 365.25)
	E = int(365.25 * D)
	G = int((C - E) / 30.6001)

	d = C - E + jd_frac - int(30.6001 * G) 
	
	if G < 13.5:
		month = G - 1
	else:
		month = G - 13
	month = int(month)
		
	if month > 2.5:
		y = D - 4716
	else:
		y = D - 4715 # year
	y = int(y)
	
	(d_frac, d) = math.modf(d)
	d = int(d)
	h = d_frac*24 # fractional part of day * 24 hours

	(h, min, s) = dec2sex(h)
	(s_frac, s) = math.modf(s)
	ms = s_frac * 60 * 1000

	s = int(s)
	ms = int(s)

	return datetime.datetime(y, month, d, h, min, s, ms)
	
def mjd2datetime(mjd):
	""" 
	Converts a Modified Julian Date to a Python datetime object. The resulting time is in UTC.
	
	Parameters
	----------
	mjd : float (any numeric type)
		A Modified Julian Date
	
	Returns
	-------
	< > : datetime.datetime
		A Python datetime.datetime object calculated using an algorithm from the book
		"Practical Astronomy with your Calculator" by Peter Duffet-Smith (Page 8)
	
	Examples
	--------
	>>> mjd2datetime(55580.90429)
	2011-01-19 21:42:10.000010
	
	"""
	jd = mjd2jd(mjd)
	return jd2datetime(jd)

def datetime2jd(datetimeObj):
	""" 
	Converts a Python datetime object to a Julian Date.
	
	Parameters
	----------
	datetimeObj : datetime.datetime
		A Python datetime.datetime object calculated using an algorithm from the book
		"Practical Astronomy with your Calculator" by Peter Duffet-Smith (Page 7)
	
	Returns
	-------
	< > : float 
		A Julian Date
	
	Examples
	--------
	>>> datetime2jd(datetimeObject)
	2455581.40429
	
	"""
	A = ymd2jd(datetimeObj.year, datetimeObj.month, datetimeObj.day)
	B = datetime2decHours(datetimeObj.time()) / 24.0
	return A + B

def datetime2mjd(datetimeObj):
	""" 
	Converts a Python datetime object to a Modified Julian Date.
	
	Parameters
	----------
	datetimeObj : datetime.datetime
		A Python datetime.datetime object calculated using an algorithm from the book
		"Practical Astronomy with your Calculator" by Peter Duffet-Smith (Page 7)
	
	Returns
	-------
	< > : float
		A Modified Julian Date
	
	Examples
	--------
	>>> datetime2mjd(datetimeObject)
	55580.90429
	
	"""
	jd = datetime2jd(datetimeObj)
	return jd2mjd(jd)

def gmst2lst(longitude, hour, minute=None, second=None, lonDirection='W', lonUnits='DEGREES'):
	""" 
	Converts Greenwich Mean Sidereal Time to Local Sidereal Time. 
	
	Parameters
	----------
	longitude : float (any numeric type)
		The longitude of the site to calculate the Local Sidereal Time. Defaults are
		Longitude WEST and units DEGREES, but these can be changed with the optional 
		parameters lonDirection and lonUnits.
	hour : int (or float)
		If an integer, the function will expect a minute and second. If a float, it
		will ignore minute and second and convert from decimal hours to hh:mm:ss.
	minute : int
		Ignored if hour is a float.
	second : int (any numeric type, to include microseconds)
		Ignored if hour is a float.
	lonDirection : string
		Default is longitude WEST, 'W', but you can specify EAST by passing 'E'.
	lonUnits : string
		Default units are 'DEGREES', but this can be switched to radians by passing
		'RADIANS' in this parameter.
	
	Returns
	-------
	hour : int
		The hour of the calculated LST
	minute : int
		The minutes of the calculated LST
	second: float
		The seconds of the calculated LST
	
	Examples
	--------
	>>> gmst2lst(70.3425, hour=14, minute=26, second=18)
	(9, 44, 55.80000000000126)
	>>> gmst2lst(5.055477, hour=14.4383333333333333, longitudeDirection='E', longitudeUnits='RADIANS')
	(9, 44, 55.79892611013463)
	
	"""
	if minute != None and second != None:
		hours = sex2dec(hour, minute, second)
	elif minute == None and second == None:
		hours = hour
	else:
		raise AssertionError('minute and second must either be both set, or both unset.')
	
	if lonUnits.upper() == 'DEGREES':
		lonTime = longitude / 15.0
	elif lonUnits.upper() == 'RADIANS':
		lonTime = longitude * 180.0 / math.pi / 15.0
	
	if lonDirection.upper() == 'W':
		lst = hours - lonTime
	elif longitudeDirection.upper() == 'E':
		lst = hours + lonTime
	else:
		raise AssertionError('longitudeDirection must be W or E')
	
	lst = lst % 24.0
	
	return dec2sex(lst)
	
def lst2gmst(longitude, hour, minute=None, second=None, lonDirection='W', lonUnits='DEGREES'):
	""" 
	Converts Local Sidereal Time to Greenwich Mean Sidereal Time.
	
	Parameters
	----------
	longitude : float (any numeric type)
		The longitude of the site to calculate the Local Sidereal Time. Defaults are
		Longitude WEST and units DEGREES, but these can be changed with the optional 
		parameters lonDirection and lonUnits.
	hour : int (or float)
		If an integer, the function will expect a minute and second. If a float, it
		will ignore minute and second and convert from decimal hours to hh:mm:ss.
	minute : int
		Ignored if hour is a float.
	second : int (any numeric type, to include microseconds)
		Ignored if hour is a float.
	lonDirection : string
		Default is longitude WEST, 'W', but you can specify EAST by passing 'E'.
	lonUnits : string
		Default units are 'DEGREES', but this can be switched to radians by passing
		'RADIANS' in this parameter.
	
	Returns
	-------
	hour : int
		The hour of the calculated GMST
	minute : int
		The minutes of the calculated GMST
	second: float
		The seconds of the calculated GMST
	
	Examples
	--------
	>>> lst2gmst(70.3425, hour=14, minute=26, second=18)
	(19, 7, 40.20000000000607)
	>>> lst2gmst(5.055477, hour=14.4383333333333333, longitudeDirection='E', longitudeUnits='RADIANS')
	(19, 7, 40.20107388985991)
	
	"""
	if minute != None and second != None:
		hours = sex2dec(hour, minute, second)
	elif minute == None and second == None:
		hours = hour
	else:
		raise AssertionError('minute and second must either be both set, or both unset.')
	
	if lonUnits.upper() == 'DEGREES':
		lonTime = longitude / 15.0
	elif lonUnits.upper() == 'RADIANS':
		lonTime = longitude * 180.0 / math.pi / 15.0
	
	if lonDirection.upper() == 'W':
		gmst = hours + lonTime
	elif longitudeDirection.upper() == 'E':
		gmst = hours - lonTime
	else:
		raise AssertionError('longitudeDirection must be W or E')
	
	gmst = gmst % 24.0
	
	return dec2sex(gmst)

def mjdLon2lst(mjd, longitude, lonDirection='W', lonUnits='DEGREES'):
	""" 
	Converts an MJD and a Longitude to Local Sidereal Time
	
	Parameters
	----------
	mjd : float (any numeric type)
		A Modified Julian Date for which to compute GMST.
	longitude : float (any numeric type)
		The longitude of the site to calculate the Local Sidereal Time. Defaults are
		Longitude WEST and units DEGREES, but these can be changed with the optional 
		parameters lonDirection and lonUnits.
	lonDirection : string
		Default is longitude WEST, 'W', but you can specify EAST by passing 'E'.
	lonUnits : string
		Default units are 'DEGREES', but this can be switched to radians by passing
		'RADIANS' in this parameter.
	
	Returns
	-------
	hour : int
		The hour of the calculated LST
	minute : int
		The minutes of the calculated LST
	second: float
		The seconds of the calculated LST
	
	"""
	dt = mjd2datetime(mjd)
	gmst = utcDatetime2gmst(dt)
	return gmst2lst(longitude, hour=gmst.hour, minute=gmst.minute, second=(gmst.second+gmst.microsecond/10**6), lonDirection=lonDirection, lonUnits=lonUnits)

#
# Coordinate Transformations:
def raLST2ha(ra, lst, raUnits='DEGREES'):
	""" 
	Converts a Right Ascension to an Hour Angle, given the Local Sidereal Time.
	
	Parameters
	----------
	ra : float, string, Angle
		A right ascension. If a float or int, default units are DEGREES but 
		can be changed to HOURS using raUnits. If a float, must be in 
		'HH:MM:SS.SS' format
	lst : float, string, Angle
		A Local Sidereal Time. If float or int, default units are HOURS. If a
		string, must be in 'HH:MM:SS.SS' format, and raUnits is ignored.
	raUnits : string
		Can be either HOURS or DEGREES, defaults to HOURS.
	
	Returns
	-------
	HA : Angle
		Returns the Hour Angle as an Angle object, which acts like a float in
		units degrees unless you want to do conversions.
	
	Notes
	-----
	For more information on the Angle object, see apwlib.astro.geometry
	
	Examples
	--------
	>>> ha = raLST2ha(70.3425, "14:26:18", raUnits='DEGREES')
	>>> print ha
	146.23250000 degrees
	>>> print ha.hours
	9.74883333333
	
	"""
	if not isinstance(ra, geometry.Angle):
		_ra = geometry.Angle()
		if isinstance(ra, "".__class__):
			_ra.hms = ra
		else:
			if raUnits == 'DEGREES':
				_ra.deg = ra
			elif raUnits == 'HOURS':
				_ra.hours = ra
			else:
				raise AssertionError('Invalid raUnits, must be HOURS or DEGREES')
	else:
		_ra = ra
	
	if not isinstance(lst, geometry.Angle):
		_lst = geometry.Angle()
		if isinstance(lst, "".__class__):
			_lst.hms = lst
		else:
			_lst.hours = lst
	else:
		_lst = lst
	
	HA = geometry.Angle()
	HA.hours = _lst.hours - _ra.hours
	
	return HA

def haLST2ra(HA, lst, haUnits='DEGREES'):
	""" 
	Converts an Hour Angle and Local Sidereal Time to Right Ascension.
	
	Parameters
	----------
	HA: float, string, Angle
		An Hour Angle. If float or int, default units are DEGREES, but can be changed
		to HOURS with haUnits. If a string, must be in 'HH:MM:SS.SS' format, 
		and raUnits is ignored.
	lst : float, string, Angle
		A Local Sidereal Time. If float or int, default units are HOURS. If a
		string, must be in 'HH:MM:SS.SS' format, and raUnits is ignored.
	haUnits : string
		Can be either HOURS or DEGREES, defaults to HOURS.
	
	Returns
	-------
	ra : Angle
		Returns the Right Ascension as an Angle object, which acts like a float in
		units degrees unless you want to do conversions.
	Notes
	-----
	For more information on the Angle object, see apwlib.astro.geometry
	
	Examples
	--------
	>>> raLST2ha(70.3425, hour=14, minute=26, second=18, raUnits='DEGREES')
	9.74883333333
	>>> raLST2ha(5.055477, hour=14.4383333333333333)
	9.74883333333
	
	"""
	
	if not isinstance(ha, geometry.Angle):
		_ha = geometry.Angle()
		if isinstance(ha, "".__class__):
			_ha.hms = ha
		else:
			if haUnits == 'DEGREES':
				_ha.deg = ha
			elif raUnits == 'HOURS':
				_ha.hours = ha
			else:
				raise AssertionError('Invalid haUnits, must be HOURS or DEGREES')
	else:
		_ha = ha
	
	if not isinstance(lst, geometry.Angle):
		_lst = geometry.Angle()
		if isinstance(lst, "".__class__):
			_lst.hms = lst
		else:
			_lst.hours = lst
	else:
		_lst = lst
	
	RA = geometry.Angle()
	RA.hours = _lst.hours - _ha.hours
	
	return RA

# TODO: 
def raUTC2ha():
	pass
	
# TODO: 
def haUTC2ra():
	pass

# TODO: 
def ra2transitTime():
	pass
	
def eclipticLatLon2RADec(lat, lon, latLonUnits='DEGREES', raDecUnits='DEGREES'):
	""" 
	Converts an Ecliptic Latitude and Ecliptic Longitude to a Right Ascension and 
	Declination.
	
	Parameters
	----------
	lat : float (any numeric type)
		An ecliptic latitude, default units DEGREES
	long : float (any numeric type)
		An ecliptic longitude, default units DEGREES	
	latLonUnits : string
		Can be either HOURS or DEGREES, defaults to DEGREES.
	raDecUnits : string
		Can be either HOURS or DEGREES, defaults to DEGREES.
		
	Returns
	-------
	ra : float
		A right ascension, default units DEGREES
	dec : float
		A declination, default units DEGREES
		
	Examples
	--------
	>>> eclipticLatLon2RADec(70.3425, -11.4552)
	(143.72252629028003, 19.535734683739964)
	
	"""	
	if latLonUnits.upper() == 'HOURS':
		lambdaa = lat * 15.0
		beta = lon * 15.0
	elif latLonUnits.upper() == 'DEGREES':
		lambdaa = lat
		beta = lon
	else:
		raise AssertionError('latLonUnits must be either HOURS or DEGREES')
	
	# Python works in Radians...
	lambdaaRad = lambdaa * math.pi / 180.0
	betaRad = beta * math.pi / 180.0
	eeRad = ee * math.pi / 180.0
	
	deltaRad = math.asin(math.sin(betaRad)*math.cos(eeRad) + math.cos(betaRad)*math.sin(eeRad)*math.sin(lambdaaRad))
	
	y = math.sin(lambdaaRad)*math.cos(eeRad) - math.tan(betaRad)*math.sin(eeRad)
	x = math.cos(lambdaaRad)
	
	alphaPrime = math.atan(y/x) * 180.0 / math.pi
	
	if y > 0 and x > 0:
		alpha = alphaPrime % 90.0
	elif y > 0 and x < 0:
		alpha = (alphaPrime % 90.0) + 90.0
	elif y < 0 and x < 0:
		alpha = (alphaPrime % 90.0) + 180.0
	elif y < 0 and x > 0:
		alpha = (alphaPrime % 90.0) + 270.0
	else:
		alpha = alphaPrime
		
	if raDecUnits.upper() == 'HOURS':
		ra = alpha / 15.0
		dec = deltaRad * 180.0 / math.pi / 15.0
	elif raDecUnits.upper() == 'DEGREES':
		ra = alpha
		dec = deltaRad * 180.0 / math.pi
	else:
		raise AssertionError('raDecUnits must be either HOURS or DEGREES')
	
	return ra, dec
	
def raDec2EclipticLatLon(ra, dec, latLonUnits='DEGREES', raDecUnits='DEGREES'):
	""" 
	Converts a Right Ascension and Declination to an Ecliptic Latitude and 
	Ecliptic Longitude.
	
	Parameters
	----------
	ra : float (any numeric type)
		A right ascension, default units DEGREES
	dec : float (any numeric type)
		A declination, default units DEGREES	
	latLonUnits : string
		Can be either HOURS or DEGREES, defaults to DEGREES.
	raDecUnits : string
		Can be either HOURS or DEGREES, defaults to DEGREES.
		
	Returns
	-------
	lat : float
		An ecliptic latitude, default units DEGREES
	lon : float
		An ecliptic longitude, default units DEGREES
		
	Examples
	--------
	>>> 
	
	
	"""
	if raDecUnits.upper() == 'HOURS':
		ra = ra * 15.0
		dec = dec * 15.0
	elif raDecUnits.upper() == 'DEGREES':
		ra = ra
		dec = dec
	else:
		raise AssertionError('raDecUnits must be either HOURS or DEGREES')
	
	# Python works in Radians...
	raRad = ra * math.pi / 180.0
	decRad = dec * math.pi / 180.0
	eeRad = ee * math.pi / 180.0
	
	betaRad = math.asin(math.sin(decRad)*math.cos(eeRad) - math.cos(decRad)*math.sin(eeRad)*math.sin(raRad))
	
	y = math.sin(raRad)*math.cos(eeRad) + math.tan(decRad)*math.sin(eeRad)
	x = math.cos(raRad)
	
	lambdaPrime = math.atan(y/x) * 180.0 / math.pi
	
	if y > 0 and x > 0:
		lambdaa = lambdaPrime % 90.0
	elif y > 0 and x < 0:
		lambdaa = (lambdaPrime % 90.0) + 90.0
	elif y < 0 and x < 0:
		lambdaa = (lambdaPrime % 90.0) + 180.0
	elif y < 0 and x > 0:
		lambdaa = (lambdaPrime % 90.0) + 270.0
	else:
		lambdaa = lambdaPrime
		
	if latLonUnits.upper() == 'HOURS':
		lat = lambdaa / 15.0
		lon = betaRad * 180.0 / math.pi * 15.0
	elif latLonUnits.upper() == 'DEGREES':
		lat = lambdaa
		lon = betaRad * 180.0 / math.pi
	else:
		raise AssertionError('latLonUnits must be either HOURS or DEGREES')
	
	return lat, lon

if __name__ == '__main__':
	import logging
	import argparse
	
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-v','--verbose', help='Be chatty.', action='store_true', dest='verbose')
	args = parser.parse_args()
	
	if args.verbose:
		logging.basicConfig(level=logging.DEBUG)
	else:
		logging.basicConfig(level=logging.INFO)
	
	# Color the output
	okColor = '\033[92m'
	errorColor = '\033[91m'
	ENDC = '\033[0m'
	
	# ----------------------------
	# The same day in MJD / JD / SDSSJD:
	mjd = 55555
	jd = 2455555.5
	sdssjd = 55555.3
	year = 2010
	month = 12
	day = 25
	hour = 0
	minute = 0
	second = 0
	theDatetime = datetime.datetime(year, month, day, hour, minute, second)
	
	# Unit Tests:
	if mjd2jd(mjd) != jd: logging.error("mjd2jd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("mjd2jd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("mjd2jd() = %s, jd = %s" % (mjd2jd(mjd), jd))
	
	if jd2mjd(jd) != mjd: logging.error("jd2mjd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("jd2mjd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("jd2mjd() = %s, mjd = %s" % (jd2mjd(jd), mjd))
	
	if mjd2sdssjd(mjd) != sdssjd: logging.error("mjd2sdssjd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("mjd2sdssjd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("mjd2sdssjd() = %s, sdssjd = %s" % (mjd2sdssjd(mjd), sdssjd))
	
	if sdssjd2mjd(sdssjd) != mjd: logging.error("sdssjd2mjd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("sdssjd2mjd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("sdssjd2mjd() = %s, mjd = %s" % (sdssjd2mjd(sdssjd), mjd))
	
	if str(jd2datetime(jd)) != str(theDatetime): logging.error("sdssjd2mjd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("sdssjd2mjd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("jd2datetime() = %s, theDatetime = %s" % (jd2datetime(jd), theDatetime))
	
	if ymd2jd(year, month, day) != jd: logging.error("ymd2jd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("ymd2jd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("ymd2jd() = %s, jd = %s" % (ymd2jd(year, month, day), jd))
	
	if datetime2jd(theDatetime) != jd: logging.error("datetime2jd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("datetime2jd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("datetime2jd() = %s, jd = %s" % (datetime2jd(theDatetime), jd))
	
	if datetime2mjd(theDatetime) != mjd: logging.error("datetime2mjd() [" + errorColor + "ERROR" + ENDC + "]")
	else: logging.info("datetime2mjd() [" + okColor + "OK" + ENDC + "]")
	logging.debug("datetime2mjd() = %s, mjd = %s" % (datetime2mjd(theDatetime), mjd))

	