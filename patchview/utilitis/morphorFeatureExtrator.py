import numpy as np
import pandas as pd
import patchview.neurom as morphor_nm

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