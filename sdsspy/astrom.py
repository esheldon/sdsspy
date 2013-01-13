"""
Defines the class Astrom for converting between pixels and equatorial
coordinates.
"""
import numpy
from numpy import where, array, cos, sin, arctan2, \
        arcsin, zeros, deg2rad, rad2deg, abs

class Astrom(object):
    """
    Conversions between pixel coordinates and equatorial coordinates.

    The intermediate step is the "great circle" coordinates mu,nu.

    dependencies
    ------------
    numpy
    You need the asTrans file for the run of interest.

    For the inverse transforms from ra,dec to pixels you need
    scipy
    """
    def __init__(self, run, camcol, filter):
        from . import util
        self._run=run
        self._camcol=camcol
        self._fchar=util.FILTERCHAR[filter]
        self._fnum=util.FILTERNUM[filter]

        self.load()


    def pix2eq(self, field, row, col, rmi=0.3):
        """
        convert from pixel coordinates (row,col) to equatorial coordinates
        (ra,dec)

        parameters
        ----------
        field: integer
            The SDSS field.
        row,col: scalar or arrays
            The pixel coords

        outputs
        -------
        ra,dec:
            The equatorial coords in degrees
        """
 
        mu,nu=self.pix2munu(field,row,col,rmi=rmi)
        ra,dec=self.gc2eq(mu,nu)
        return ra,dec

    def eq2pix(self, field, ra, dec, rmi=0.3):
        """
        convert from equatorial coordinates (ra,dec) to pixel coordinates
        (row,col) 

        parameters
        ----------
        field: integer
            The SDSS field.
        ra,dec:
            The equatorial coords in degrees

        outputs
        -------
        row,col: scalar or arrays
            The pixel coords
        """
 
        mu,nu = self.eq2gc(ra, dec)
        row,col = self.munu2pix(field, mu, nu, rmi=rmi)

        return row,col

    def eq2gc(self, ra, dec):
        """
        convert from equatorial (ra,dec) to SDSS great circle coordinates
        (mu,nu)

        parameters
        ----------
        field: integer
            The SDSS field.
        ra,dec: scalar or arrays
            The equatorial coordinates in degrees

        outputs
        -------
        mu,nu:
            The SDSS great circle coords in degrees
        """
        ra,dec,are_scalar=self._get_array_args(ra,dec,"ra","dec")

        rarad=deg2rad(ra)
        decrad=deg2rad(dec)
 
        node=self._node_rad
        incl=self._incl_rad

        cosdec=cos(decrad)
        x1 = cos(rarad - node)*cosdec
        y1 = sin(rarad - node)*cosdec
        z1 = sin(decrad)
        x2 = x1
        y2 = y1*cos(incl) + z1*sin(incl)
        z2 =-y1*sin(incl) + z1*cos(incl)

        mu = arctan2(y2, x2) + node
        nu = arcsin(z2)

        rad2deg(mu,mu)
        rad2deg(nu,nu)

        self._atbound2(nu, mu)

        if are_scalar:
            mu=mu[0]
            nu=nu[0]

        return mu,nu

    def gc2eq(self, mu, nu):
        """
        convert from SDSS great circle coordinates (mu,nu) to equatorial
        (ra,dec)

        parameters
        ----------
        field: integer
            The SDSS field.
        mu,nu: scalar or arrays
            The great circle coordinates in degrees

        outputs
        -------
        ra,dec:
            Equatorial coordinates in degrees
        """
 
        mu,nu,are_scalar=self._get_array_args(mu,nu,"mu","nu")

        murad=deg2rad(mu)
        nurad=deg2rad(nu)

        node=self._node_rad
        incl=self._incl_rad
        
        cosnu=cos(nurad)
        x2 = cos(murad-node)*cosnu
        y2 = sin(murad-node)*cosnu
        z2 = sin(nurad)
        y1 = y2*cos(incl) - z2*sin(incl)
        z1 = y2*sin(incl) + z2*cos(incl)

        ra = arctan2(y1, x2) + node
        dec = arcsin(z1)

        rad2deg(ra,ra)
        rad2deg(dec,dec)

        self._atbound2(dec, ra)

        if are_scalar:
            ra=ra[0]
            dec=dec[0]

        return ra, dec

    def munu2pix(self, field, mu, nu, rmi=0.3):
        """
        Solve for the row,col that are roots of the equation

            row,col -> mu,nu

        This is because we only have the forward transform

        parameters
        ----------
        field: integer
            SDSS field number
        mu,nu: 
            SDSS great circle coordinates in degrees

        outputs
        -------
        row,col:
            pixel values
        """

        import scipy.optimize
        mu,nu,are_scalar=self._get_array_args(mu,nu,"mu","nu")

        rmi=array(rmi,dtype='f8', ndmin=1, copy=False)
        if rmi.size > 1 and rmi.size !=mu.size:
            raise ValueError("rmi should be scalar or same size as inputs")

        trans=self._astrans
        w,=where(trans['field']==field)
        if w.size == 0:
            raise ValueError("field not found: %s" % field)

        # pack away this data for the optimizer
        fd={'a':trans['a'][w],
            'b':trans['b'][w],
            'c':trans['c'][w],
            'd':trans['d'][w],
            'e':trans['e'][w],
            'f':trans['f'][w],
            'drow0':trans['drow0'][w],
            'drow1' : trans['drow1'][w],
            'drow2' : trans['drow2'][w],
            'drow3' : trans['drow3'][w],

            'dcol0' : trans['dcol0'][w],
            'dcol1' : trans['dcol1'][w],
            'dcol2' : trans['dcol2'][w],
            'dcol3' : trans['dcol3'][w],

            'csrow' : trans['csrow'][w],
            'cscol' : trans['cscol'][w],
            'ccrow' : trans['ccrow'][w],
            'cccol' : trans['cccol'][w],
            'ricut' : trans['ricut'][w]}

        self._field_data=fd

        det = fd['b']*fd['f'] - fd['c']*fd['e']
        mudiff = mu - fd['a']
        nudiff = nu - fd['d']
        row_guess = ( mudiff*fd['f'] - fd['c']*nudiff )/det
        col_guess = ( fd['b']*nudiff - mudiff*fd['e'] )/det
 
        row=zeros(mu.size,dtype='f8')
        col=zeros(mu.size,dtype='f8')
        for i in xrange(mu.size):
            if rmi.size==1:
                self._tmp_rmi=rmi[0]
            else:
                self._tmp_rmi=rmi[i]

            self._tmp_munu=array([mu[i],nu[i]])

            rowcol_guess=array([row_guess[i], col_guess[i]])

            rowcol = scipy.optimize.fsolve(self._pix2munu_for_fit, rowcol_guess)
            row[i] = rowcol[0]
            col[i] = rowcol[1]

        if are_scalar:
            row=row[0]
            col=col[0]

        return row,col


    def pix2munu(self, field, row, col, rmi=0.3):
        """
        convert from row-column to SDSS great circle coordinates (mu,nu)

        parameters
        ----------
        field: integer
            The SDSS field.
        row,col: scalar or arrays
            The row,column in the field.

        outputs
        -------
        mu,nu:
            SDSS great circle coords in degrees
        """
       
        row,col,are_scalar=self._get_array_args(row,col,"row","col")

        trans=self._astrans
        w,=where(trans['field']==field)
        if w.size == 0:
            raise ValueError("field not found: %s" % field)

        a = trans['a'][w]
        b = trans['b'][w]
        c = trans['c'][w]
        d = trans['d'][w]
        e = trans['e'][w]
        f = trans['f'][w]

        drow0 = trans['drow0'][w]
        drow1 = trans['drow1'][w]
        drow2 = trans['drow2'][w]
        drow3 = trans['drow3'][w]

        dcol0 = trans['dcol0'][w]
        dcol1 = trans['dcol1'][w]
        dcol2 = trans['dcol2'][w]
        dcol3 = trans['dcol3'][w]

        csrow = trans['csrow'][w]
        cscol = trans['cscol'][w]
        ccrow = trans['ccrow'][w]
        cccol = trans['cccol'][w]
        ricut = trans['ricut'][w]
      
        rowm = 0.5 + row+drow0+drow1*col+drow2*(col**2)+drow3*(col**3)+csrow*rmi
        colm = 0.5 + col+dcol0+dcol1*col+dcol2*(col**2)+dcol3*(col**3)+cscol*rmi

        mu = a + b * rowm + c * colm
        nu = d + e * rowm + f * colm

        if are_scalar:
            mu=mu[0]
            nu=nu[0]

        return mu,nu

    def load(self):
        import esutil as eu
        from . import files
        

        fname=files.filename('asTrans', run=self._run, 
                             camcol=self._camcol)
        h=eu.io.read_header(fname,ext=0)
        self._header=h
        self._node=h['node']
        self._incl=h['incl']
        self._node_rad=deg2rad(self._node)
        self._incl_rad=deg2rad(self._incl)

        camcols=array([int(c) for c in h['camcols'].split()])
        filters=array(h['filters'].split())

        wcam,=where(camcols==self._camcol)
        if wcam.size==0:
            raise  ValueError("bad camcol: %s" % self._camcol)

        wfilter,=where(filters==self._fchar)
        if wfilter.size==0:
            raise  ValueError("bad filter: %s" % self._fchar)
        wcam,wfilter=wcam[0],wfilter[0]

        ext = wcam*len(filters) + (wfilter + 1)

        trans,hthis=eu.io.read(fname,ext=ext,lower=True,header=True)

        thisfilt=hthis['filter'].strip()
        thiscam=hthis['camcol']
        if self._fchar != thisfilt:
            raise ValueError("read wrong filter: %s instead of "
                             "%s" % (thisfilt,self._fchar))
        if self._camcol != thiscam:
            raise ValueError("read wrong camcol: %s instead of "
                             "%s" % (thiscam,self._camcol))

        self._astrans=trans


    def _pix2munu_for_fit(self, rowcol):
        """
        This is a version for the fitter only
        """
        row=rowcol[0]
        col=rowcol[1]
        fd=self._field_data

        rowm = (0.5 + row +
                fd['drow0']
                +fd['drow1']*col
                +fd['drow2']*(col**2)
                +fd['drow3']*(col**3)
                +fd['csrow']*self._tmp_rmi)

        colm = (0.5 + col
                +fd['dcol0']
                +fd['dcol1']*col
                +fd['dcol2']*(col**2)
                +fd['dcol3']*(col**3)
                +fd['cscol']*self._tmp_rmi)

        diff=zeros(2,dtype='f8')
        diff[0] = fd['a'] + fd['b'] * rowm + fd['c'] * colm
        diff[1] = fd['d'] + fd['e'] * rowm + fd['f'] * colm

        diff[0] = abs(diff[0]-self._tmp_munu[0])
        diff[1] = abs(diff[1]-self._tmp_munu[1])

        return diff

    def _atbound(self, angle, minval, maxval):

        w,=where(angle < minval)
        if w.size > 0:
            angle[w] = angle[w] + 360.0
            w,=where(angle < minval)

        w,=where(angle >= maxval)
        while w.size != 0:
            angle[w] = angle[w] - 360.0
            w,=where(angle >= maxval)


    def _atbound2(self, theta, phi):

        self._atbound( theta, -180.0, 180.0 )
        w, = where( abs(theta) > 90.)
        if w.size > 0:
            theta[w] = 180. - theta[w]
            phi[w] = phi[w] + 180.

        self._atbound( theta, -180.0, 180.0 )
        self._atbound( phi, 0.0, 360.0 )

        w,=where( abs(theta) == 90. )
        if w.size > 0:
            phi[w] = 0.0

    def _get_array_args(self, x1, x2, name1, name2):
        is_scalar1=_is_scalar(x1)
        is_scalar2=_is_scalar(x2)

        if is_scalar1 != is_scalar2:
            mess="%s and %s must both be scalar or array"
            mess = mess % (name1,name2)
            raise ValueError(mess)

        x1=array(x1,copy=False,dtype='f8',ndmin=1)
        x2=array(x2,copy=False,dtype='f8',ndmin=1)

        if x1.size != x2.size:
            mess="%s and %s must both be same length.  got %s and %s"
            mess = mess % (name1,name2,x1.size,x2.size)
            raise ValueError(mess)

        return x1,x2,is_scalar1

def _is_scalar(obj):
    try:
        l=len(obj)
        return False
    except:
        return True
