import numpy as np
from scipy import stats
import pandas as pd
import neurom as morphor_nm
from neurom import morphmath as mm
from neurom.core.types import (tree_type_checker, NeuriteType)
from neurom import features
from scipy.spatial import ConvexHull, Delaunay
from scipy.spatial.distance import cdist
from neurom.core.morphology import Section

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

def extractMorhporFeatures_simple(n, df_summary=None):
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
    df_summary["Number of neurite"] = [number_of_neurites]
    df_summary["Number of sections"] = [number_of_sections]
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
    segment_length = np.array(morphor_nm.get('segment_lengths', n, neurite_type=neurite_type))
    if len(segment_midpoints)==0:
        return [], [], [], []     
    n_density = np.array(segment_midpoints)[:, 'xyz'.index(axis)]
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
    x_pts = np.array(segment_midpoints)[:, 'xyz'.index('x')]
    y_pts = np.array(segment_midpoints)[:, 'xyz'.index('y')]
    if useFullRange == False:      
        px_min, px_max = np.min(x_pts), np.max(x_pts)     
        py_min, py_max = np.min(y_pts), np.max(y_pts) 
        binX = int((px_max - px_min)/step_size)    
        binY = int((py_max - py_min)/step_size) 
    else:
        segment_midpoints = morphor_nm.get('segment_midpoints', n, neurite_type=NeuriteType.all)
        segment_length_ = morphor_nm.get('segment_lengths', n, neurite_type=NeuriteType.all)
        x_pts_ = np.array(segment_midpoints)[:, 'xyz'.index('x')]     
        y_pts_ = np.array(segment_midpoints)[:, 'xyz'.index('y')]
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

def getPerimeter_Area(ps):
    ps2  = ps.copy()
    ps2 = ps2[:,:2] ## make sure we get 2D array
    try:
        c_hull = ConvexHull(ps2)
        area, perimeter =  c_hull.volume,c_hull.area
    except:
        print('Hull alg. failed! return None for area and perimeter for the soma')
        area, perimeter =  np.nan, np.nan
    return area, perimeter

def getSomaCircularityIndex(area, perimeter):
    ''' Calculate roundness of a soma contour 
    return a float value between (0,1.0). 1 being a perfect circle
    '''
    if area is not np.nan and perimeter is not np.nan:
        return  np.round((4 * np.pi * area) / (perimeter * perimeter), 2)
    else:
        return np.nan

def getShapeFactors(ps_):
    ''' 
    Giving soma points,return a bunch of shape factors in one go
    circularity, max_diameter, shape_factor as defined in neurom.morphmath
    '''
    ps = ps_[:,:2].copy()
    try:
        hull = ConvexHull(ps)
        cirIndex = 4.0 * np.pi * hull.volume / hull.area**2 #circularity
        hull_points = ps[hull.vertices]
        max_pairwise_distance = np.max(cdist(hull_points, hull_points))
        shapefactor= hull.volume / max_pairwise_distance**2
    except:
        cirIndex, max_pairwise_distance, shapefactor  = np.nan, np.nan, np.nan
    try:
        aspRatio = mm.aspect_ratio(ps)
    except:
        aspRatio = np.nan
    return cirIndex, max_pairwise_distance, shapefactor,aspRatio

def sec_len(n,sec):
    """Return the length of a section."""
    return mm.section_length(sec.points)

def total_neurite_length(n):
    '''Total neurite length (sections)'''
    return sum(sec_len(n, s) for s in morphor_nm.iter_sections(n))

def total_neurite_volume(n):
    '''Total neurite volume'''
    return sum(mm.segment_volume(s) for s in morphor_nm.iter_segments(n))

def total_neurite_area(n):
    '''Total neurite area'''
    return sum(mm.segment_area(s) for s in morphor_nm.iter_segments(n))

def total_bifurcation_points(n):
    '''Number of bifurcation points'''
    return sum(1 for _ in morphor_nm.iter_sections(n,
                                          iterator_type=Section.ibifurcation_point))
def max_branch_order(n):
    '''Maximum branch order'''
    return  max(features.section.branch_order(s) for s in morphor_nm.iter_sections(n))

def min_max_trunk_angle(n, neurteType = NeuriteType.basal_dendrite):
    angles = morphor_nm.features.morphology.trunk_angles(n, neurite_type=neurteType, consecutive_only=False)
    if len(angles)>1:
        # print(n.name, angles)
        nmax = []
        nmin = []
        for x in angles:
            if 0 in x:
                x.remove(0)
            nmax.append(max(x))
            nmin.append(min(x))
        try:
            ad_min = min(nmin)
            ad_max = max(nmax)
        except:
            ad_min = np.nan
            ad_max = np.nan    
    else:
        ad_min = np.nan
        ad_max = np.nan    
    return ad_min, ad_max, angles

def dendrites_dispersion_index(ad_min, ad_max, angles):
    '''dendrites_dispersion_index'''
    ## calculate dendrite directions dispersion
    ## ad_max: maximum angle between any two dendrites, in (0, pi)
    ## ad_min: minimum angle between any two dendrites, in (0, pi)
    ## dispersion index: (ad_max-ad_min) / (ad_max+ad_min) (0, 1). 0 being total overlap.
    if ad_max!=np.nan and ad_min!=np.nan:
        if ad_max==ad_min:
            ad_max = np.pi ## two dendrites only
            return 1-(ad_max-ad_min) / (ad_max+ad_min)
    else:
        return np.nan

def extractMorphFeatures(n, df_summary=None):
    ''' return a dictionary contains useful morphor features
    n: NeuroM Neuron object
    '''
    if df_summary is None:
        df_summary = {}
    maxDia, soma_center, soma_radius, soma_avgRadius = getSomaStats(n)
    df_summary["Neuron id"] = n.name
    df_summary["center X"] = soma_center[0]
    df_summary["center Y"] = soma_center[1]
    df_summary["center Z"] = soma_center[2]
    cellArea, cellPerimeter = getPerimeter_Area(n.soma.points)
    cirIndex, max_pairwise_distance, shapefactor, asratio = getShapeFactors(n.soma.points)
    df_summary["soma average radius"] = soma_avgRadius
    df_summary["soma maximal radius"] = np.max(soma_radius)
    df_summary["soma minimal radius"] = np.min(soma_radius)
    df_summary['soma max_pairwise_dist'] = max_pairwise_distance
    df_summary["soma perimeter"] = cellPerimeter
    df_summary["soma area"] = cellArea
    df_summary['soma circularity index'] = cirIndex
    df_summary['soma shape factor'] = shapefactor
    df_summary['soma aspect_ratio'] = asratio
    df_summary['cell max_radial_dist'] = features.morphology.max_radial_distance(n)
    df_summary['total number of neurites'] = len(n.neurites)
    if len(n.neurites) > 0:
        for neurite in n.neurites:
            nsections = features.morphology.number_of_sections_per_neurite(n, neurite.type)
            neuriteTrunkLength = morphor_nm.features.morphology.trunk_section_lengths(n, neurite_type=neurite.type)
            neutriteName = str(neurite.type).split('.')[-1]
            df_summary[neutriteName+ ' Nseg'] = len(nsections)
            tortuositys = [morphor_nm.features.section.section_tortuosity(s) for s in morphor_nm.iter_sections(n, 
            neurite_filter=tree_type_checker(neurite.type))]
            df_summary[neutriteName+ '_avg_tortuosity'] = np.mean(tortuositys)
            df_summary[neutriteName+ '_max_tortuosity'] = np.max(tortuositys)
            if len(nsections) > 1 :
                for i in range(len(nsections)):
                    df_summary[neutriteName+ ' seg'+str(i+1)+' Nsec'] = nsections[i]
                    df_summary[neutriteName+ ' seg'+str(i+1)+' trunkLen'] = neuriteTrunkLength[i] # stem length
                
            else:
                df_summary[neutriteName+' Nsec'] = nsections[0]
                df_summary[neutriteName+' trunkLen'] = neuriteTrunkLength[0] # stem length
        # neurite features
        neurite_funcs = [total_neurite_length,total_neurite_volume,total_neurite_area,\
            total_neurite_area,total_bifurcation_points,
            max_branch_order]
        for f in neurite_funcs:
            df_summary[f.__doc__] = f(n)

    admin, admax, angles = min_max_trunk_angle(n)
    df_summary['trunk angle min'] = admin
    df_summary['trunk angle max'] = admax
    df_summary['trunk angle dispersion index'] = dendrites_dispersion_index(admin, admax, angles)
    return df_summary

if __name__ == "__main__":
    pass