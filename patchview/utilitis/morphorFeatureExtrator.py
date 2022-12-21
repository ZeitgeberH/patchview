import numpy as np
from scipy import stats
import pandas as pd
import neurom as morphor_nm
from neurom.core.types import NeuriteType
from neurom import features
def getSomaStats(n):
    ''' Basic statistics for soma
    '''
    ps = n.soma.points[:, :3]
    maxDia = minmaxDist(ps)
    ps_center = np.nanmean(ps, axis=0)
    ps_radius = np.sqrt(np.nansum((ps - ps_center) ** 2, axis=1))
    ps_avgRadius = np.nanmean(ps_radius)
    return maxDia, ps_center, ps_radius, ps_avgRadius

def minmaxDist(ps):
    """calculate minimal and maximal diameter
    ps: pointset from a neuron object. Type: 2-D numpy array
    """
    nPoints = len(ps)
    # minDia = 1000
    maxDia = 0
    for j in range(nPoints - 1):
        for k in range(j + 1, nPoints):
            diameter = np.sqrt(np.sum((ps[j] - ps[k]) ** 2))
            # if diameter < minDia:
            #     minDia = diameter
            if diameter > maxDia:
                maxDia = diameter
    return maxDia

def extractMorhporFeatures(n, df_summary=None):
    ''' return a dataframe (formated as record) contains useful features
    '''
    number_of_neurites = morphor_nm.get("number_of_neurites", n)
    # Extract the total number of sections
    number_of_sections = morphor_nm.get("number_of_sections", n)
    # Extract the number of sections per neurite
    number_of_sections_per_neurite = morphor_nm.get(
        "number_of_sections_per_neurite", n
    )
    number_of_sections_per_neurite = morphor_nm.get(
    "number_of_sections_per_neurite", n)

    if df_summary is None:
        df_summary = {}
    maxDia, soma_center, soma_radius, soma_avgRadius = getSomaStats(n)
    df_summary["Neuron id"] = [n.name]
    df_summary["center X"] = [soma_center[0]]
    df_summary["center Y"] = [soma_center[1]]
    df_summary["center Z"] = [soma_center[2]]
    df_summary["average radius"] = [soma_avgRadius]
    df_summary["maximal radius"] = [np.max(soma_radius)]
    df_summary["minimal radius"] = [np.min(soma_radius)]
    df_summary["maximal diameter"] = [maxDia]
    df_summary["Number of neurite"] = [number_of_neurites[0]]
    df_summary["Number of sections"] = [number_of_sections[0]]
    for i, neurite in enumerate(n.neurites):
        df_summary[str(neurite.type)] = [number_of_sections_per_neurite[i]]
    # df = pd.DataFrame(df_summary)
    df = pd.DataFrame.from_dict(df_summary, orient="index").transpose()
    df = np.array(
        df.to_records(index=False)
    )  ## format dataframe for using in QtTable widget
    return df

def sholl_analysis(n, step_size=10,neurite_type=NeuriteType.all):
    ''' sholl analysis
    n: NeuroM neuron object with neurites
    step_size: step size for analysis
    output: frequency and bins (in um)
    '''
    freq = morphor_nm.get("sholl_frequency", n, neurite_type=neurite_type, step_size=step_size)
    bins =  list(n.soma.radius + np.arange(len(freq))*step_size)
    return freq, bins

def sholl_single_axis(n, step_size=1, axis='x', neurite_type=NeuriteType.all):
    ''' projection to single axis.
        Return density (um per length step). distance is center at soma center.
    '''
    segment_midpoints = morphor_nm.get('segment_midpoints', n, neurite_type=neurite_type)
    segment_length = morphor_nm.get('segment_lengths', n, neurite_type=neurite_type)
    if len(segment_midpoints)==0:
        return [], [], [], []      
    n_density = np.array(segment_midpoints[:, 'xyz'.index(axis)]) 
    maxDia, soma_center, soma_radius, soma_avgRadius = getSomaStats(n)
    ps_min, ps_max = np.min(n_density), np.max(n_density)
    if axis == 'x':
        c = soma_center[0]
    else:
        c = soma_center[1]
    nBins = int((ps_max - ps_min)/step_size)
    seg_counts, bin_edges = np.histogram(n_density, bins = nBins, density=False) ## number of segments
    seg_len_idx = np.digitize(n_density, bin_edges) ## get bin index for each segment
    density = np.zeros(len(bin_edges))
    for n in range(len(bin_edges)): ## get total segment length within each bin
        density[n] = np.sum(segment_length[seg_len_idx==n])    
    return density, bin_edges, c, n_density

def sholl_2D_density(n,  step_size=2, neurite_type=NeuriteType.all, maxNorm=False, useFullRange=False):
    if type(neurite_type) is not list:
        segment_midpoints = morphor_nm.get('segment_midpoints', n, neurite_type=neurite_type)
        segment_length = morphor_nm.get('segment_lengths', n, neurite_type=neurite_type)
    else:
        segment_midpoints = []
        segment_length = []
        for nt in neurite_type:
            segment_midpoints.extend(morphor_nm.get('segment_midpoints', n, neurite_type=nt))
            segment_length.extend(morphor_nm.get('segment_lengths', n, neurite_type=nt))
        segment_midpoints = np.array(segment_midpoints)
        segment_length = np.array(segment_length)

    if len(segment_midpoints)==0:
        return [], [], [], [], []
    x_pts = np.array(segment_midpoints[:, 'xyz'.index('x')])     
    y_pts = np.array(segment_midpoints[:, 'xyz'.index('y')])
    if useFullRange == False:      
        px_min, px_max = np.min(x_pts), np.max(x_pts)     
        py_min, py_max = np.min(y_pts), np.max(y_pts) 
        binX = int((px_max - px_min)/step_size)    
        binY = int((py_max - py_min)/step_size) 
    else:
        segment_midpoints = morphor_nm.get('segment_midpoints', n, neurite_type=NeuriteType.all)
        segment_length_ = morphor_nm.get('segment_lengths', n, neurite_type=NeuriteType.all)
        x_pts_ = np.array(segment_midpoints[:, 'xyz'.index('x')])     
        y_pts_ = np.array(segment_midpoints[:, 'xyz'.index('y')])
        px_min, px_max = np.min(x_pts_), np.max(x_pts_)     
        py_min, py_max = np.min(y_pts_), np.max(y_pts_) 
        nbinX = int((px_max - px_min)/step_size)    
        nbinY = int((py_max - py_min)/step_size) 
        ret = stats.binned_statistic_2d(x_pts_, y_pts_, segment_length_, 'sum', bins=[nbinX, nbinY], expand_binnumbers=True)
        binX = ret.x_edge
        binY = ret.y_edge
    ret = stats.binned_statistic_2d(x_pts, y_pts, segment_length, 'sum', bins=[binX, binY], expand_binnumbers=True)
    maxDia, soma_center, soma_radius, soma_avgRadius = getSomaStats(n)
    if maxNorm:
        dataStatistcs = ret.statistic/np.max(ret.statistic)
    else:
        dataStatistcs = ret.statistic
    return dataStatistcs, ret.x_edge, ret.y_edge, soma_center[0], soma_center[1]

def sholl_polar(n,  step_size=1, pho_step=5, angle_step=np.pi/16,  neurite_type=NeuriteType.all):
    dataStats, x_edge, y_edge, centerX, centerY = sholl_2D_density(n, step_size=step_size, neurite_type=neurite_type, maxNorm=False)
    polar_data = dataStats.reshape((-1,1))
    polar_coords_P = []
    polar_coords_A = []
    c = 0
    x_edge -=centerX
    y_edge -=centerY
    for idx, x in enumerate(x_edge[1:]):
        for idy, y in enumerate(y_edge[1:]):
            if polar_data[c] > 0:
                polar_coords_P.append(np.sqrt(x**2 + y**2))
                polar_coords_A.append(np.arctan2(y,x))
            c += 1
    rbins = np.linspace(0, np.max(polar_coords_P), int(np.max(polar_coords_P)//pho_step), endpoint=False)
    abins = np.linspace(-np.pi,np.pi, int(2*np.pi//angle_step), endpoint=True)          
    hist, _, _ = np.histogram2d(np.array(polar_coords_A), np.array(polar_coords_P), bins=(abins, rbins))
    A, R = np.meshgrid(abins, rbins)
    return hist, A, R


def sholl_2D_density_(n, step_size=2, neurite_type=NeuriteType.all, bins=None, maxNorm=False):
    ''' density in x-y plane.
        Return density (um per length step). distance is center at soma center.
        TODO: USE scipy.stats.binned_statistic_2d to get density
    '''
    horiz, bin_h, centerH, segH = sholl_single_axis(n, step_size=step_size, axis='x', neurite_type=neurite_type)
    vert, bin_v, centerV, segV = sholl_single_axis(n, step_size=step_size, axis='y', neurite_type=neurite_type)
    # ndiff =  bin_h.shape[0] - bin_v.shape[0]
    # if ndiff>0:
    #     bin_v = padding_(ndiff, bin_v, padVal=None)
    #     vert = padding_(ndiff, vert, padVal=0)
    # else:
    #     bin_h = padding_(-ndiff, bin_h, padVal=None)
    #     horiz = padding_(-ndiff, horiz, padVal=0)
    if bins is None:
        bins = [bin_h, bin_v]
    else:
        bins = [bins, bins]
    if maxNorm:
        d2d, h_edges, v_edges = np.histogram2d(segH, segV, bins=bins, density=False)
        d2d = d2d/np.max(d2d)
    else:
        d2d, h_edges, v_edges = np.histogram2d(segH, segV, bins=bins, density=True)

    return d2d, h_edges, v_edges, centerH, centerV

if __name__ == "__main__":
    pass