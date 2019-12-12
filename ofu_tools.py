import numpy as np

# from the flarestack software:
def ofu_ang_distance(evt1,evt2):
    return ofu_angular_distance(evt1['ra'],evt1['dec'],evt2['ra'],evt2['dec'])

def ofu_angular_distance(lon1, lat1, lon2, lat2):
    """
    calculate the angular distince along the great circle
    on the surface of a shpere between the points
    (`lon1`,`lat1`) and (`lon2`,`lat2`)
    This function Works for equatorial coordinates
    with right ascension as longitude and declination
    as latitude. This function uses the Vincenty formula
    for calculating the distance.
    Parameters
    ----------
    lon1 : array_like
      longitude of first point in radians
    lat1 : array_like
      latitude of the first point in radians
    lon2 : array_like
      longitude of second point in radians
    lat2 : array_like
      latitude of the second point in radians
    """
    c1 = np.cos(lat1)
    c2 = np.cos(lat2)
    s1 = np.sin(lat1)
    s2 = np.sin(lat2)
    sd = np.sin(lon2 - lon1)
    cd = np.cos(lon2 - lon1)

    return np.arctan2(
        np.hypot(c2 * sd, c1 * s2 - s1 * c2 * cd),
        s1 * s2 + c1 * c2 * cd
    )


def ofu_doublet_average(start,end,theta_a=0.9):
    # theta_a : telescope FoV term

    # Compute the average direction of doublet, weighted with 1/uncertainty**2
    ra = np.column_stack((start['ra'],end['ra']))
    dec = np.column_stack(( start['dec'],end['dec']))
    weight = np.column_stack((1./start['angErr']**2,1./start['angErr']**2))

    x = np.average(weight * np.cos(ra) * np.cos(dec),axis=1)
    y = np.average(weight * np.sin(ra) * np.cos(dec),axis=1)
    z = np.average(weight * np.sin(dec),axis=1)
    r_average = np.sqrt(x**2. + y**2. + z**2.)
    ra_average = np.arctan2(y, x)
    # bring ra into [0,2*pi]
    ra_average = np.where(ra_average < 0.,2.*np.pi - abs(ra_average),ra_average)
    dec_average = np.arcsin(z/(r_average))
    # directional error is: dirError = totalSigma, totalSigma=sqrt(1/sumOfWeights), weights are 1/sigma^2
    ang_error = np.sqrt(1/np.sum(weight,axis=1))

    # calculate TS:
    ang_d = ofu_angular_distance(start['ra'],start['dec'],end['ra'],end['dec'])
    time_d = end['time']-start['time'] #in MJD -> in TS: convert 100s to MJD
    sigma_q_2 = start['angErr']**2+end['angErr']**2
    # test statistic according according to OFU paper:
    TS = ang_d**2/sigma_q_2 + 2*np.log(2*np.pi*sigma_q_2) - 2*np.log(1-np.exp(-theta_a**2/(2*ang_error))) + 2*np.log(time_d/(100/86400.))
    return np.rec.fromarrays( (ra_average,dec_average,ang_error,start['time'],end['time'],TS), names=('ra_av', 'dec_av', 'sigma_av', 't0', 't1', 'ts') )
    #return np.ndarray( [(ra_average,dec_average,ang_error,end['time'],TS)],dtype=[('ra_av', 'f8'), ('dec_av', 'f8'), ('sigma_av', 'f8'), ('end_time', 'f8'), ('ts', 'f8')] )


# general purpose
def ofu_ComputeAverageDirection(data):
    #Compute the average direction and uncertainty of all events in data, weighted with 1/uncertainty**2
    ra = data['ra']
    dec = data['dec']
    weight = 1./data['sigma']**2

    x = np.average(weight * np.cos(ra) * np.cos(dec))
    y = np.average(weight * np.sin(ra) * np.cos(dec))
    z = np.average(weight * np.sin(dec))
    r_average = np.sqrt(x**2. + y**2. + z**2.)
    ra_average = np.arctan2(y, x)
    if ra_average < 0.:
        ra_average = 2.*np.pi - abs(ra_average)
    dec_average = np.arcsin(z/(r_average))
    # directional error is: dirError = totalSigma, totalSigma=sqrt(1/sumOfWeights), weights are 1/sigma^2
    ang_error = np.sqrt(1/sum(weight))

    return np.array( [(ra_average,dec_average,ang_error)],dtype=[('ra_av', 'f8'), ('dec_av', 'f8'), ('sigma_av', 'f8')] )
