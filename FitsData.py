import os

import fitsio

class FitsData:
    """
    Collects FITS file data specified at file_path so that the FITS file 
    may be closed. Also Provides methods to organize collected FITS data.

    Based on fitsio.FITS class (uses cfitsio)

    Parameters:
        file_path: The path to the FITS file to be opened

    Attributes:
        file_path: file_path provided at instantiation
        file_size: size of file at file_path
        data: List of HDU units
        data_bytype: List of dict data tables by table extension name
    """

    def __init__(self, file_path):

        # File related attributes
        self.file_path = file_path
        self.file_size = str(os.stat(file_path).st_size)

        def _collect_data(self, file_path):
            """Collects data from FITS file"""

            self.data = []
            with fitsio.FITS(file_path) as file_handle:
                for hdu in file_handle:
                    header = hdu.read_header()

                    if 'NAXIS2' in hdu.read_header() \
                        and hdu.read_header()['NAXIS2'] is not 0:
                        table = hdu.read()   
                    else:
                        table = []

                    self.data.append({'header': header, 'table': table})

        def _organize_data(self):
            """Organizes data by extension table name"""

            self.data_bytype = {}
            for hdu in self.data:
                if 'EXTNAME' in hdu['header']:
                    extname = hdu['header']['EXTNAME'].strip()

                    # Create entry
                    if extname not in self.data_bytype:
                        self.data_bytype[extname] = []

                    self.data_bytype[extname].append(hdu['table'])

        _collect_data(self, file_path)
        _organize_data(self)

    def __getitem__(self, ind):
        """Returns HDU at index"""

        return self.data[ind]






















            





