import netCDF4 as nc
import numpy as np
import os

class Database(object):
    """
    Dict-like data structure for accessing NetCDF variables contained in multiple
    files. It allows variable lookup using variable names alone, which is convenient
    as long as distinct variables have unique lines even when comparing across different
    files. It provides no guarantees about relationship between variables in different
    files---for example, variables in different files that share a dimension name are not
    guaranteed to share the same coordinate value, or even to have the same size along
    that dimension.

    xarray (xarray.pydata.org) provides a much more complete set of tools for working
    with multi-file netCDF datasets. This class is not intended to duplicate that
    functionality---its role is to reduce boilerplate in code that needs the low-level 
    control provided by the netCDF4 library.
    """

    def __init__(self, files_to_vars):
        """
        Create a variable database from a dict-like mapping from file
        paths to variable names. This is used to create an inverse
        mapping from variable names to open files.

        Parameters
        ----------
        files_to_vars: dict-like (str -> list of str)
            Mapping from file paths to variable names. Files
            paths must point to netCDF files that contain
            each variable in the corresponding list; if not,
            an exception will be raised
        """

        self.table = dict()
        
        # Try block catches exceptions in order to close
        # already-opened files before re-raising them
        try:
            for path, var_list in files_to_vars.items():
                ncfile = nc.Dataset(path, 'r')
                for var in var_list:
                    if not var in ncfile.variables.keys():
                        ncfile.close()
                        raise ValueError(
                            'NetCDF file at %s does not contain variable %s' % (path, var))
                    if var in self.table.keys():
                        ncfile.close()
                        raise ValueError(
                            'Attempting to add duplicate entry %s to database' % var)
                    self.table[var] = ncfile

        except err:
            self.close()
            raise err

    def close(self):
        """
        Close all file handles stored in the database
        """
        for var, ncfile in self.table.items():
            if ncfile.isopen():
                ncfile.close()

    def lookup(self, var):
        """
        Look up and return a variable from a NetCDF file in the database

        self[var] can be used as a shortcut for self.lookup(var)

        Parameters
        ----------
        var: str
            Variable key

        Returns
        -------
        netCDF4.Variable
            NetCDF variable corresponding to the provided key
        """
        return self.table[var].variables[var]

    def files(self):
        """
        Return a list of netCDF file paths
        """
        files = []
        for var, ncfile in self.table.items():
            if not ncfile.filepath() in files:
                files.append(ncfile.filepath())
        return files

    def coordinates(self, var):
        """
        Return a dict mapping the names of a variable's dimensions to the underlying
        netCDF coordinates (or None, if the dimension is unlabeled).
        """
        dset = self.table[var]
        dimensions = dset.variables[var].dimensions
        coords = dict(zip(
            dimensions, 
            (dset.variables[dim] if dim in dset.variables else None for dim in dimensions)
        ))
        return coords

    def __getitem__(self, var):
        return self.lookup(var)

    def __str__(self):
        summary = "<class 'quicknc.Database'>"
        for var, ncfile in self.table.items():
            summary += '\n%s: %s' % (var, self[var].__str__())
        return summary

def create(path, file_metadata):
    """
    Create a writeable netCDF4 Dataset based on provided metadata. This function
    primarily exists to reduce boilerplate when creating new netCDF files; see the
    netCDF documentation for details about file metadata.

    Parameters
    ----------
    path: str
        Location of the newly-created file. Trying to create a file that already exists 
        will overwrite the existing file.
    
    file_metadata: dict-like
        File metadata. Only three keys---'vars', 'coords', and 'attrs'---are
        used by this function. The corresponding items should contain the following

        vars: dict-like
            Mapping from variable name to datatype, list of coordinate names and 
            dict of attributes. Coordinate names must be correspond to keys in the 
            'coords' item.

        coords: dict-like
            Mapping from coordinate names to datatype, coordinate length and dict of attributes. 
            Coordinate lengths can be integers or None; the latter corresponds to an 
            unlimited dimension. 
            
            A one-dimensional array can also be provided in place of the coordinate datatype and length, 
            in which case the array is used to define the coordinate datatype and size and to fill the
            corresponding variable.

        attrs: dict-like
            Global attributes for the netCDF file

    Returns
    -------
    netCDF4.Dataset
        Dataset for the new netCDF file, empty except (possibly) for provided coordinate arrays.
    """

    # Create file
    head, tail = os.path.split(os.path.abspath(path))
    if not os.path.exists(head):
        os.makedirs(head, exist_ok=True)
    dset = nc.Dataset(path, mode='w', clobber=True)

    # Create coordinates 
    # This is most of the work
    for name, data in file_metadata['coords'].items():
        
        if len(data) == 3:
            dtype, size, attrs = data
        elif len(data) == 2:
            arr, attrs = data
            if arr.ndim != 1:
                raise ValueError('Coordinate array for %s is not one-dimensional' % name)
            dtype = arr.dtype
            size = arr.size
        else:
            raise TypeError('Could not parse metadata for coordinate %s' % name)

        dset.createDimension(name, size=size)
        var = dset.createVariable(name, dtype, dimensions=(name,))
        if len(data) == 2:
            var[:] = arr[:]

        for key, value in attrs.items():
            setattr(var, key, value)

    # Create variables
    for name, (dtype, coords, attrs) in file_metadata['vars'].items():    
        var = dset.createVariable(name, dtype, dimensions=coords)
        for key, value in attrs.items():
            setattr(var, key, value)

    # Set global attributes
    for key, value in file_metadata['attrs'].items():
        setattr(dset, key, value)

    return dset
        
            

        
        
        
    