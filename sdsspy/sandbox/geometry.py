################################################################################
#
# geometry.py - A module containing various Python objects for performing
#				Geometric calculations.
#
# Copyright (C) 2011 Adrian Price-Whelan
# 
################################################################################

"""
todo:
"""

__author__ = 'Adrian Price-Whelan <adrn@nyu.edu>'

# Standard library dependencies (e.g. sys, os)
from math import radians, degrees, pi, fabs

class Angle(object):
	@property
	def deg(self):			
		return self._deg
	
	@deg.setter
	def deg(self, angle):
		if isinstance(angle, Angle):
			angle = angle.deg
		self._deg = float(angle)
		self._rad = float(radians(self._deg))
		self._hours = float(self._deg / 15.0)
		try:
			self.hour, self.minute, self.second = dec2sex(self._hours)
		except NameError:
			from apwlib.astro.convert import dec2sex
			self.hour, self.minute, self.second = dec2sex(self._hours)
		self._hms = "%02d:%02d:%02.2f" % (self.hour, self.minute, self.second)
	
	@property
	def rad(self):
		return self._rad
	
	@rad.setter
	def rad(self, angle):
		if isinstance(angle, Angle):
			angle = angle.rad
		self._rad = float(angle)
		self._deg = float(degrees(self._rad))
		self._hours = float(self._deg / 15.0)
		try:
			self.hour, self.minute, self.second = dec2sex(self._hours)
		except NameError:
			from apwlib.astro.convert import dec2sex
			self.hour, self.minute, self.second = dec2sex(self._hours)
		self._hms = "%02d:%02d:%02.2f" % (self.hour, self.minute, self.second)
	
	@property
	def hours(self):
		return self._hours
	
	@hours.setter
	def hours(self, angle):
		if isinstance(angle, Angle):
			angle = angle.hours
		self._hours = float(angle)
		self._deg = float(self._hours * 15.0)
		self._rad = float(radians(self._deg))
		try:
			self.hour, self.minute, self.second = dec2sex(self._hours)
		except NameError:
			from apwlib.astro.convert import dec2sex
			self.hour, self.minute, self.second = dec2sex(self._hours)
		self._hms = "%02d:%02d:%02.2f" % (self.hour, self.minute, self.second)
	
	@property
	def hms(self):
		return self._hms
	
	@hms.setter
	def hms(self, angle):
		if isinstance(angle, Angle):
			angle = angle.hms
		self._hms = angle
		try:
			self._hours = string2hours(angle)
		except NameError:
			from apwlib.astro.convert import string2hours
			self._hours = string2hours(angle)
		self._deg = float(self._hours * 15.0)
		self._rad = float(radians(self._deg))
		try:
			self.hour, self.minute, self.second = dec2sex(self._hours)
		except NameError:
			from apwlib.astro.convert import dec2sex
			self.hour, self.minute, self.second = dec2sex(self._hours)
	
	def __float__(self):
		return self._deg
	
	def __int__(self):
		return int(self._deg)
	
	def __str__(self):
		return "%.8f%s" % (self._deg, " degrees")
	
	def __add__(self, other):
		return self._deg + other
	
	def __radd__(self, other):
		return other + self._deg
	
	def __sub__(self, other):
		return self._deg - other
	
	def __rsub__(self, other):
		return other - self._deg
	
	def __mul__(self, other):
		return self._deg * other
	
	def __rmul__(self, other):
		return other * self._deg
	
	def __div__(self, other):
		return self._deg / other
	
	def __rdiv__(self, other):
		return other / self._deg
	
	def __radians__(self):
		return self._rad
		
if __name__ == '__main__':
	ra = Angle()
	ra.hours = 9.16123344
	print ra
	print ra + 56.0
	print 56.0 + ra
	print ra - 56.0
	print 56.0 - ra
	print ra / 56.0
	print 56.0 / ra
	print ra * 56.0
	print 56.0 * ra
	