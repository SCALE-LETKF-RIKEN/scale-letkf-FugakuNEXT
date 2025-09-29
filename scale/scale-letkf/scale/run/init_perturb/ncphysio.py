import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

dimlist_default = [['time'], ['z'], ['y'], ['x']]


def ncphys_read_dimlen(rootgrp, dimlist=dimlist_default):
    """
    Read lengths of dimensions.

    Parameters
    ----------
    rootgrp : netcdf4-python Dataset instance
        The input NetCDF file.
    dimlist : array of array, optional
        List of dimensions in the NetCDF file. Default: `dimlist_default`

    Returns
    -------
    dimlen : dictionary
        Lengths of dimensions in the NetCDF file.
    """
    dimlen = {}
    for idiml in dimlist:
        for idim in idiml:
            if idim in rootgrp.dimensions:
                dimlen[idim] = len(rootgrp.dimensions[idim])
            else:
                dimlen[idim] = None
    return dimlen


def ncphys_read(rootgrp, varname, dimlist=dimlist_default, time=None, it=None):
    """
    Read a variable from a NetCDF file.

    Can choose to read a single time or all times.
    The dimensions of the variable are reordered based on the 'dimlist' setting.

    Parameters
    ----------
    rootgrp : netcdf4-python Dataset instance
        The input NetCDF file.
    varname : string
        The variable name.
    dimlist : array of array, optional
        List of dimensions in the NetCDF file. Default: `dimlist_default`
    time : number, optional
        Target time in physical time unit.
    it : int, optional
        Target time index. Defalut to 0 if both `time` and `it` are not given.

    Returns
    -------
    vardim : dictionary
        Dimensions of the return variable data.
    vardata : ndarray or masked_array
        Variable data in a ndarray or masked_array (if the variable has the 
        `_FillValue` attribute).
    """
    # Check if the variable exists in the NetCDF file
    if varname not in rootgrp.variables:
        raise IOError("Variable '" + varname + "' does not exist.")

    # Check the variable dimensions and the order of the variable dimensions
    # in the input NetCDF file
    vardim = rootgrp.variables[varname].dimensions
    dimord = [None] * len(dimlist)
    dimname = [None] * len(dimlist)
    for i, idiml in enumerate(dimlist):
        for j, idim in enumerate(vardim):
            if idim in idiml:
                if dimord[i] is None:
                    dimord[i] = j
                    dimname[i] = idim
                else:
                    raise ValueError('Duplicated dimensions.')
    if len([i for i in dimord if i is not None]) != len(vardim):
        raise ValueError("Variable has addtional dimensions not listed in 'dimlist'")

    # Read the variable data into a ndarray or masked_array
    if dimord[0] is None or time == 'all' or it == 'all':
        vardata = rootgrp.variables[varname][:]
    else:
        if time is not None:
            found = np.where(rootgrp.variables[dimname[0]][:] == time)[0]
            if len(found) == 0:
                raise ValueError("Cannot find 'time' = " + str(time))
            else:
                it0 = found[0]
        elif it is None:
            it0 = 0
        elif it < 0 or it >= len(rootgrp.dimensions['time']):
            raise ValueError("Cannot find 'it' = " + str(it))
        else:
            it0 = it
        if dimord[0] == 0:
            vardata = rootgrp.variables[varname][it0]
        else:
            vardata = np.take(rootgrp.variables[varname][:], it0, axis=dimord[0])
    if '_FillValue' in rootgrp.variables[varname].ncattrs() and type(vardata) != ma.MaskedArray:
        vardata = ma.array(vardata, fill_value=rootgrp.variables[varname].getncattr('_FillValue'))

    # Reorder the axes of the return variable data if necessary
    if time == 'all' or it == 'all':
        vardim_out = [i for i in dimname if i is not None]
        transaxes = [i for i in dimord if i is not None]
    else:
        vardim_out = [i for i in dimname[1:] if i is not None]
        transaxes = []
        for iord in dimord[1:]:
            if iord is not None:
                if dimord[0] is not None and iord > dimord[0]:
                    transaxes.append(iord - 1)
                else:
                    transaxes.append(iord)
    if transaxes != list(range(len(transaxes))):
        vardata = np.transpose(vardata, axes=transaxes)

    return vardim_out, vardata


def ncphys_write(rootgrp, varname, vardim, vardata, dimlist=dimlist_default, time=None, it=None):
    """
    Write a variable to a NetCDF file.

    Can choose to write a single time or all times.

    Parameters
    ----------
    rootgrp : netcdf4-python Dataset instance
        The input NetCDF file.
    varname : string
        The variable name.
    vardim : dictionary
        Dimensions of the input variable data.
    vardata : ndarray or masked_array
        Variable data to be written to the files.
    dimlist : array of array, optional
        List of dimensions in the NetCDF file. Default: `dimlist_default`
    time : number, optional
        Target time in physical time unit.
    it : int, optional
        Target time index. Defalut to 0 if both `time` and `it` are not given.
    """
    # Check if the variable exists in the NetCDF file
    if varname in rootgrp.variables:
        vardim_in = rootgrp.variables[varname].dimensions
    else:
        raise IOError("Variable '" + varname + "' does not exist.")
#        ncphys_create_var(vardim)
#        vardim_in = vardim

    # Check the variable dimensions and the order of the variable dimensions
    # in the input NetCDF file, and reorder the axes of the data if necessary
    dimord = [None] * len(vardim_in)
    tdim = None
    for i, idim in enumerate(vardim_in):
        if idim in vardim:
            dimord[i] = vardim.index(idim)
            if vardata.shape[dimord[i]] != len(rootgrp.dimensions[idim]):
                raise ValueError("Variable dimensions mismatch.")
        elif idim in dimlist[0] and tdim is None and (time != 'all' and it != 'all'):
            dimord[i] = None
            tdim = i
            if time is not None:
                found = np.where(rootgrp.variables[idim][:] == time)[0]
                if len(found) == 0:
                    raise ValueError('Cannot find time = ' + str(time))
                else:
                    it0 = found[0]
            elif it is None:
                it0 = 0
            elif it < 0 or it >= len(rootgrp.dimensions['time']):
                raise ValueError("Cannot find 'it' = " + str(it))
            else:
                it0 = it
        else:
            raise ValueError("Variable dimensions mismatch.")
    transaxes = [i for i in dimord if i is not None]
    if len(transaxes) != len(vardim):
        raise ValueError("Variable dimensions mismatch.")
    if transaxes != list(range(len(transaxes))):
        vardata = np.transpose(vardata, axes=transaxes)

    # Write the variable data
    if tdim is None:
        rootgrp.variables[varname][:] = vardata
    else:
        if tdim == 0:
            rootgrp.variables[varname][it0] = vardata
        else:
            slice_obj = [slice(None)] * len(vardim_in)
            slice_obj[tdim] = it0
            rootgrp.variables[varname][slice_obj] = vardata
