#
# In FITS files, the DATASUM and CHECKSUM keywords play an important role in data integrity verification. 
# The DATASUM  keyword holds the sum of the data values.
# The CHECKSUM keyword holds an ASCII-encoded checksum of the Header Data Unit (HDU).
#
# A detailed description of these keywords and the algorithms used for calculating them can be found in 
# the FITS standard document:
# https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
#
# In the following, I provide an example of how to create or update these keywords using Python and the astropy library.
#
# It is crucial to note that any updates to the DATASUM and CHECKSUM keywords must be performed as the last 
# operation before closing the file. This is because any subsequent modification to the data in the file 
# (such as changing a keyword or altering the data itself) will render the DATASUM and CHECKSUM values 
# inconsistent with the data. In such cases, the DATASUM and CHECKSUM values must be recalculated and updated.

#!/usr/bin/python
﻿from astropy.io import fits
with fits.open('name_of_your_file.fits', mode='update') as list_hdu:
    # Loop through each HDU
    for hdu in list_hdu:
        # Add the DATASUM keyword
        hdu.add_datasum()
        # Add the CHECKSUM keyword
        hdu.add_checksum()
    # Save the changes to the file
    hdu_list.flush()
