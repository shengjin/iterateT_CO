"""
    readImage() - Reads RADMC3D image(s)
                - Module needed for readimage_fft.py
"""

try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '

from copy import deepcopy



##########################################*********************************************
class radmc3dImage():
    """
    RADMC3D image class

    ATTRIBUTES:
    -----------
        image       - The image as calculated by radmc3d (the values are intensities in erg/s/cm^2/Hz/ster)
        imageJyppix - The image with pixel units of Jy/pixel
        x           - x coordinate of the image [cm]
        y           - y coordinate of the image [cm]
        nx          - Number of pixels in the horizontal direction
        ny          - Number of pixels in the vertical direction
        sizepix_x   - Pixel size in the horizontal direction [cm]
        sizepix_y   - Pixel size in the vertical direction [cm]
        nfreq       - Number of frequencies in the image cube
        freq        - Frequency grid in the image cube
        nwav        - Number of wavelengths in the image cube (same as nfreq)
        wav         - Wavelength grid in the image cube

    """
    def __init__(self):
        self.image       = 0
        self.imageJyppix = 0
        self.x           = 0
        self.y           = 0
        self.nx          = 0
        self.ny          = 0
        self.sizepix_x   = 0
        self.sizepix_y   = 0
        self.nfreq       = 0
        self.freq        = 0
        self.nwav        = 0
        self.wav         = 0
        self.stokes      = False
        self.psf         = {}
        self.fwhm        = []
        self.pa          = 0
        self.dpc         = 0

#####################################################
    def readImage(self, fname=None):
        """
        Function to read an image calculated by RADMC3D 
     
        INPUT:
        ------
         fname   : file name of the radmc3d output image (if omitted 'image.out' is used)
        """

        pc   = 3.0856775e18

        # Look for the image file
        if (fname==None): 
            fname = 'image.out'

        try:
            rfile = open(fname, 'r')
        except:
            print 'ERROR!'
            print 'No '+fname+' file has been found!'
            return -1
    
        dum = ''

        # Format number
        # iformat = 1 means an observer at infinity
        iformat = int(rfile.readline())
        if iformat != 1:
            print "Error: the program only support iformat = 1 case!"

        # Nr of pixels
        dum = rfile.readline()
        dum = dum.split()
        self.nx  = int(dum[0])
        self.ny  = int(dum[1])
        print self.nx, self.ny

        # Nr of frequencies
        # self.nfreq = 1 means just one image at 1 wavelength
        self.nfreq = int(rfile.readline())
        self.nwav  = self.nfreq 
        print self.nwav, self.nfreq

        # Pixel sizes
        # the size of pixels in cm (for an observer at infinity) 
        # for infinity mode, pixel-in-arcsec = pixel-in-cm/1.496e13/distance-in-parsec 
        # or radian (for local observer mode)
        dum = rfile.readline()
        dum = dum.split()
        self.sizepix_x = float(dum[0])
        self.sizepix_y = float(dum[1])
        print self.sizepix_x, self.sizepix_y

        # Wavelength of the image
        # wavelenth in micron belonging to the various images in this file
        self.wav = []
        for iwav in range(self.nwav):
            self.wav.append(float(rfile.readline()))
        self.wav = np.array(self.wav)
        self.freq = 2.99792458e10 / self.wav * 1e4
        print self.wav, self.freq

        # We have a normal total intensity image
        self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=float64)
        for iwav in range(self.nwav):
            # Blank line
            dum = rfile.readline()
            for ix in range(self.nx):
                for iy in range(self.ny):
                    self.image[ix,iy,iwav] = float(rfile.readline())
           
        rfile.close()

        # Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
        conv  = self.sizepix_x * self.sizepix_y / pc**2. * 1e23 
        self.imageJyppix = self.image * conv
        self.x = ((arange(self.nx, dtype=float64) + 0.5) - self.nx/2) * self.sizepix_x
        self.y = ((arange(self.ny, dtype=float64) + 0.5) - self.ny/2) * self.sizepix_y


# --------------------------------------------------------------------------------------------------
    def imConv(self, fwhm=None, pa=None, dpc=1.):
        """
        Function to convolve a radmc3d image with a two dimensional Gaussian psf 
    
        INPUT:
        ------
              fwhm    : A list of two numbers; the FWHM of the two dimensional psf along the two principal axes
                            The unit is assumed to be arcsec 
              pa      : Position angle of the psf ellipse (counts from North counterclockwise)
              dpc     : Distance of the source in pc
    
        OUTPUT:
        -------
              result  : same  
              'cimage': The convolved image with the psf (unit is erg/s/cm/cm/Hz/ster)
              'image' : The original unconvolved image (unit is erg/s/cm/cm/Hz/ster)
              'psf'   : Two dimensional psf
              'x'     : first coordinate axis of the psf/image
              'y'     : second coordinate axis of the psf/image
        """
# --------------------------------------------------------------------------------------------------
# Natural constants    
        au = 1.496e13
        pc = 3.0857200e+18
    
        nx = self.nx
        ny = self.ny
        dx = self.sizepix_x / au / dpc
        dy = self.sizepix_y/au/dpc
        nfreq = self.nfreq
    
    
# Calculate the Gaussian psf
        dum   = getPSF(nx=self.nx, ny=self.ny, fwhm=fwhm, pa=pa, pscale=[dx,dy])
        psf   = dum['psf'] 
        f_psf = np.fft.fft2(psf)

        if self.stokes:
            if self.nfreq==1:
                cimage = zeros([self.nx,self.ny,4], dtype=float64)
                for istokes in range(4):
                    imag = self.image[:,:,istokes]
                    f_imag  = np.fft.fft2(imag)
                    f_cimag = f_psf * f_imag
                    cimage[:,:,istokes] = abs(np.fft.ifftshift(np.fft.ifft2(f_cimag)))
            else:
                cimage = zeros([self.nx,self.ny,4,self.nfreq], dtype=float64)
                for ifreq in range(nfreq):
                    for istokes in range(4):
                        imag = self.image[:,:,istokes,ifreq]
                        f_imag  = np.fft.fft2(imag)
                        f_cimag = f_psf * f_imag
                        cimage[:,:,istokes,ifreq] = abs(np.fft.ifftshift(np.fft.ifft2(f_cimag)))

        else:
            cimage = zeros([self.nx,self.ny,self.nfreq], dtype=float64)
            for ifreq in range(nfreq):
                imag = self.image[:,:,ifreq]
                f_imag  = np.fft.fft2(imag)
                f_cimag = f_psf * f_imag
                cimage[:,:,ifreq] = abs(np.fft.ifftshift(np.fft.ifft2(f_cimag)))


        #cimage = squeeze(cimage)
  
# Return the convolved image (copy the image class and replace the image attribute to the convolved image)

        res             = deepcopy(self)
        conv            = res.sizepix_x * res.sizepix_y / (dpc*pc)**2./ (fwhm[0] * fwhm[1] * pi / (4.*log(2.)))
        res.image       = cimage * conv
        res.imageJyppix = res.image * 1e23
        res.psf         = psf
        res.fwhm        = fwhm
        res.pa          = pa
        res.dpc         = dpc


        return res





###############################
# --------------------------------------------------------------------------------------------------
def getPSF(nx=None, ny=None, fwhm=None, pa=None, pscale=None):
    """
    Function to generate a two dimensional Gaussian PSF
    
    INPUT:
    ------
          nx      : image size in the first dimension
          ny      : image size in the second dimension
          fwhm    : full width at half maximum of the psf in each dimension [fwhm_x, fwhm_y]
          pa      : position angle of the gaussian if the gaussian is not symmetric
          pscale  : pixelscale of the image, if set fwhm should be in the same unit, if not set unit of fwhm is pixels

    OUTPUT:
    -------
          result  : dictionary containing the following keys
          'psf'   : two dimensional numpy array containing the normalized psf
          'x'     : first coordinate axis of the psf
          'y'     : seonc coordinate axis of the psf
          
    """
# --------------------------------------------------------------------------------------------------

# Create the two axes

    if (pscale!=None):
        dx,dy = pscale[0], pscale[1]
    else:
        dx,dy = 1., 1.

    x = (np.arange(nx, dtype=float64) - nx/2) * dx
    y = (np.arange(ny, dtype=float64) - ny/2) * dy

# Calculate the standard deviation of the Gaussians
    sigmax = fwhm[0] / (2.0 * sqrt(2.0 * log(2.)))
    sigmay = fwhm[1] / (2.0 * sqrt(2.0 * log(2.)))
    norm   = 1./(2. * pi * sigmax * sigmay)


# Pre-compute sin and cos angles

    sin_pa = np.sin(pa/180.*pi - pi/2.)
    cos_pa = np.cos(pa/180.*pi - pi/2.)

# Define the psf
    psf = np.zeros([nx,ny], dtype=float64)
    for ix in range(nx):
        for iy in range(ny):
            xx = cos_pa * x[ix] - sin_pa * y[iy]
            yy = sin_pa * x[ix] + cos_pa * y[iy]

            psf[ix,iy] = exp(-0.5*xx*xx/sigmax/sigmax - 0.5*yy*yy/sigmay/sigmay)

    
    # Normalize the PSF 
    psf = psf / norm 

    res = {'psf':psf, 'x':x, 'y':y}

    return res

#############################################################
def readImage(fname=None):
    """
    Function to read an image calculated by RADMC3D 
     
    INPUT:
    ------
        fname   : file name of the radmc3d output image (if omitted 'image.out' is used)
    """

    dum = radmc3dImage()
    dum.readImage(fname=fname)
    return dum


