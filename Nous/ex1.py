# *- coding: UTF-8 -*

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy import ndimage

x = np.linspace(0,2*np.pi,100)
sine = np.sin(x)

im = sine * sine[...,None]
d1 = ndimage.gaussian_filter(im, sigma=5, order=1, mode='wrap')
d2 = ndimage.gaussian_filter(im, sigma=5, order=2, mode='wrap')

##a = plt.figure(1)
##plt.subplot(131)
##plt.imshow(im)
##plt.title('original')
##plt.show()
##
##b = plt.figure(2)
##plt.subplot(132)
##plt.imshow(d1)
##plt.title('first derivative')
##plt.show()

c = plt.figure(3)
c = c.add_subplot(133)
plt.title('second derivative')
plt.imshow(d2)
plt.show()
