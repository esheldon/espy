import numpy as np
import re

surveyCenterRa  =  np.double(185.0)
surveyCenterDec =   np.double(32.5)

deg2rad = np.pi/180.
rad2deg = 180./np.pi

node = surveyCenterRa - 90.
node = node * deg2rad

etaPole = surveyCenterDec * deg2rad

__doc__="""
    A set of transformation functions for use with SDSS coordinate systems.

    eq2csurvey(): Convert between equatorial and corrected survey coords.
    csurvey2eq(): Convert between corrected survey and equatorial coords.

    Don't use these unless you have to:
        eq2survey(): Convert between equatorial and uncorrected survey coords.
        survey2eq(): Convert between uncorrected survey and equatorial coords.

    Adapted from astrotools
        Erin Sheldon, NYU, 2006-03-11
    Force data type and allow selection of dtype through keyword.
        Erin Sheldon, NYU, 2007-05-23
        
"""
# utility functions
def atbound(longitude, minval, maxval):
    w, = np.where(longitude < minval)
    while w.size > 0:
        longitude[w] += 360.0
        w, = np.where( longitude < minval )

    w, = np.where(longitude > maxval)
    while w.size > 0:
        longitude[w] -= 360.0
        w, = np.where( longitude > maxval )

    return

def atbound2(theta, phi):

    atbound(theta, -180.0, 180.0)

    w, = np.where( np.abs(theta) > 90.0 )
    if w.size > 0:
        theta[w] = 180.0 - theta[w]
        phi[w] += 180.0

    atbound(theta, -180.0, 180.0)
    atbound(phi, 0.0, 360.0)

    w, = np.where( np.abs(theta) == 90.0 )
    if w.size > 0:
        phi[w] = 0.0
           

def eq2csurvey(ra_in, dec_in, dtype=np.float64):
    """
    NAME:
      eq2csurvey
    PURPOSE:
       Convert from ra, dec to the corrected clambda, ceta 
       SDSS survey coordinate system.  It is corrected so that the
       longitude eta ranges from [-180.0, 180.0] and the latitude
       lambda ranges from [-90.0,90.0].  The standard lambda/eta 
       both range from [-180.0,180.0] which doesn't make sense.
       NOTE: lambda is often referred to as longitude but this
       is incorrect since it has poles at [-90,90]

    CALLING SEQUENCE:
      import sdss_transform
      (clambda, ceta) = sdss_transform.eq2csurvey(ra, dec, dtype=float64)

    INPUTS: 
      ra: Equatorial latitude in degrees. 
      dec: Equatorial longitude in degrees. 
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is float64. See 
        numpy.typeDict for a list of possible types.
        dtype: The data type of output.  Default is numpy.float64.

    OUTPUTS: 
      clambda: Corrected Survey longitude (actually lattitude) in degrees
      ceta: Corrected Survey latitude (actually logitude) in degrees
      
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    ra = np.array(ra_in, ndmin=1, copy=True, dtype=dtype)
    dec = np.array(dec_in, ndmin=1, copy=True, dtype=dtype)

    if (ra.size != dec.size):
        raise ValueError("RA, DEC must be same size")

    # range checking
    if (ra.min() < 0.0) | (ra.max() > 360.0):
        raise ValueError('RA must we within [0,360]')
    if (dec.min() < -90.0) | (dec.max() > 90.0):
        raise ValueError('DEC must we within [-90,90]')

    ra *= deg2rad
    dec *= deg2rad
    ra -= node

    # generate x,y,z on unit sphere, clearing memory as we go
    cdec = np.cos(dec)
    
    x = np.cos(ra)*cdec
    y = np.sin(ra)*cdec

    ra = 0; cdec = 0 # mem 
    
    z = np.sin(dec)
    
    dec = 0 # mem

    # generate clambda, ceta
    # do things in place to save memory

    # clambda = -arcsin( x ) (not a copy clambda=x)
    np.arcsin(x, x); clambda=x
    clambda *= -1

    
    # ceta = arctan2( z, y ) - etaPole
    np.arctan2( z, y, z); ceta = z
    ceta -= etaPole

    clambda *= rad2deg
    ceta *= rad2deg

    atbound(ceta, -180.0, 180.0)
    
    return (clambda, ceta)

def csurvey2eq(clambda_in, ceta_in, dtype=np.float64):
    """
    NAME:
      csurvey2eq
    PURPOSE:
       Convert corrected clambda, ceta SDSS survey coordinate system t
       equatorial coords.  

    CALLING SEQUENCE:
      import sdss_transform
      (ra, dec) = sdss_transform.csurvey2eq(clambda, ceta, dtype=float64)

    INPUTS: 
      clambda: Corrected Survey longitude (actually lattitude) in degrees
      ceta: Corrected Survey latitude (actually logitude) in degrees
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is float64. See 
        numpy.typeDict for a list of possible types.

    OUTPUTS: 
      ra: Equatorial latitude in degrees. 
      dec: Equatorial longitude in degrees. 
      
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """
    
    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    clambda = np.array(clambda_in, ndmin=1, copy=True, dtype=dtype)
    ceta = np.array(ceta_in, ndmin=1, copy=True, dtype=dtype)

    # range checking
    if (clambda.min() < -90.0) | (clambda.max() > 90.0):
        raise ValueError('CLAMBDA must we within [-90,90]')
    if (ceta.min() < -180.0) | (ceta.max() > 180.0):
        raise ValueError('CETA must we within [-180,180]')

    clambda *= deg2rad
    ceta    *= deg2rad

    x = -np.sin(clambda)
    y = np.cos(ceta + etaPole)*np.cos(clambda)
    z = np.sin(ceta + etaPole)*np.cos(clambda)
    
    ra = np.arctan2( y, x ) + node
    dec = np.arcsin(z)

    ra  *= rad2deg
    dec *= rad2deg
    atbound2(dec, ra)

    return (ra,dec)


def eq2survey(ra_in, dec_in, dtype=np.float64):
    """
    NAME:
      eq2survey
    PURPOSE:
       Convert from ra, dec to the lambda, eta 
       SDSS survey coordinate system.  Note this coordinate system is
       not well defined.  Recommend you use csurvey coords.

    CALLING SEQUENCE:
      import sdss_transform
      (lambda, eta) = sdss_transform.eq2survey(ra, dec, dtype=float64)

    INPUTS: 
      ra: Equatorial latitude in degrees. 
      dec: Equatorial longitude in degrees. 
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is float64. See 
        numpy.typeDict for a list of possible types.

    OUTPUTS: 
      lambda: SDSS Survey longitude (actually lattitude) in degrees
      eta: SDSS Survey latitude (actually logitude) in degrees
      
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    # Make a copy as an array. ndmin=1 to avoid messed up scalar arrays
    ra = np.array(ra_in, ndmin=1, copy=True, dtype=dtype)
    dec = np.array(dec_in, ndmin=1, copy=True, dtype=dtype)

    if (ra.size != dec.size):
        raise ValueError("RA, DEC must be same size")

    # range checking
    if (ra.min() < 0.0) | (ra.max() > 360.0):
        raise ValueError('RA must we within [0,360]')
    if (dec.min() < -90.0) | (dec.max() > 90.0):
        raise ValueError('DEC must we within [-90,90]')

    ra *= deg2rad
    dec *= deg2rad
    ra -= node

    # generate x,y,z on unit sphere, clearing memory as we go
    cdec = np.cos(dec)
    
    x = np.cos(ra)*cdec
    y = np.sin(ra)*cdec

    ra = 0; cdec = 0 # mem 
    
    z = np.sin(dec)
    
    dec = 0 # mem

    # generate lam, eta
    # do things in place to save memory

    # lam = -arcsin( x ) (not a copy lam=x)
    np.arcsin(x, x); lam=x
    lam *= -1

    
    # eta = arctan2( z, y ) - etaPole
    np.arctan2( z, y, z); eta = z
    eta -= etaPole

    lam *= rad2deg
    eta *= rad2deg

    atbound2(lam, eta)
    atbound(eta, -180.0, 180.0)

    w, = np.where( eta > (90.0 - surveyCenterDec) )
    if w.size > 0:
        eta[w] -= 180.0
        lam[w]  = 180.0 - lam[w]

    atbound(lam, -180.0, 180.0)

    return (lam, eta)


def survey2eq(ra, dec, dtype=np.float64):
    """
    NAME:
      survey2eq
    PURPOSE:
       Convert clambda, ceta SDSS survey coordinate system to
       equatorial coords.  

    CALLING SEQUENCE:
      import sdss_transform
      (ra, dec) = sdss_transform.survey2eq(lam, eta, dtype=float64)

    INPUTS: 
      lambda: Survey longitude (actually lattitude) in degrees
      eta:    Survey latitude (actually logitude) in degrees
    OPTIONAL INPUTS:
        dtype: The data type of output.  Default is float64. See 
        numpy.typeDict for a list of possible types.
      
    OUTPUTS: 
      ra: Equatorial latitude in degrees. 
      dec: Equatorial longitude in degrees. 
      
    REVISION HISTORY:
      Written: 11-March-2006  Converted from IDL program.
    """

    return csurvey2eq(ra,dec, dtype=dtype)
