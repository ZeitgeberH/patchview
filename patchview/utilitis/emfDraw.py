from morphio import SomaType
import neurom as morphor_nm
from neurom import viewer as morphor_viewer
from neurom.io.multiSomas import MultiSoma
from neurom.core.morphology import Morphology
from neurom.core.types import NeuriteType
import numpy as np
import pyemf3 as pyemf
from neurom.core.morphology import iter_neurites, iter_sections, iter_segments
from neurom.core.types import tree_type_checker
from neurom.morphmath import segment_radius
from neurom.core.dataformat import COLS
from scipy.spatial import ConvexHull
TREE_COLOR = {NeuriteType.basal_dendrite: (0x0f, 0x23, 0xbd), #'blue',
              NeuriteType.apical_dendrite: (0xcc,0x06,0xf9), #'purple',
              NeuriteType.axon: (0xf9,0x06, 0x18), #'red',
              NeuriteType.soma: (0x00,0x00,0x00), #black,
              NeuriteType.undefined: (0x81,0x91,0x59),## light green
              NeuriteType.custom6: (0x00,0x10,0x00),
              NeuriteType.custom6: (0x00,0x00,0x00),
              NeuriteType.custom7: (0x00,0x00,0x00),
              NeuriteType.custom8:(0x00,0x00,0x00),
              NeuriteType.custom9: (0x00,0x00,0x00),
              NeuriteType.custom10: (0x00,0x00,0x00)}
def _plane2col(plane):
    """Take a string like 'xy', and return the indices from COLS.*."""
    planes = ('xy', 'yx', 'xz', 'zx', 'yz', 'zy')
    assert plane in planes, 'No such plane found! Please select one of: ' + str(planes)
    return (getattr(COLS, plane[0].capitalize()),
            getattr(COLS, plane[1].capitalize()), )

def _get_linewidth(tree, linewidth, diameter_scale, DiameterScaling=2.0):
    """Calculate the desired linewidth based on tree contents.

    If diameter_scale exists, it is used to scale the diameter of each of the segments
    in the tree
    If diameter_scale is None, the linewidth is used.
    """
    if diameter_scale is not None and tree:
        linewidth = [DiameterScaling * segment_radius(s) * diameter_scale
                     for s in iter_segments(tree)]
    return linewidth

def _get_color(treecolor, tree_type):
    """If treecolor set, it's returned, otherwise tree_type is used to return set colors."""
    if treecolor is not None:
        return treecolor
    return TREE_COLOR.get(tree_type, (0x00,0x00,0x00))

def emf_soma(emf, soma,plane='xy', scalingFactor = 100):
    plane0, plane1 = _plane2col(plane)
    points = np.vstack([soma.points[:,plane0].ravel(),
                        soma.points[:,plane1].ravel()])
    points = points.T
    if points.shape[0] >0:
        try:
            hull = ConvexHull(points)        
            pt_list = (np.array(points[hull.vertices])*scalingFactor).astype('int32')
            penWidth = 1
        except Exception as e:
            pt_list = np.array(points*scalingFactor).astype('int32')
            penWidth = int(soma.radius*scalingFactor)
            print('Error in ConvexHull', penWidth)
        pt_list[:,1] = -pt_list[:,1]
        somaColor = TREE_COLOR.get(NeuriteType.soma, (0x00,0x00,0x00))
        pen = emf.CreatePen(pyemf.PS_SOLID, penWidth,somaColor)
        emf.SelectObject(pen)
        brush=emf.CreateSolidBrush(somaColor)
        emf.SelectObject(brush)
        emf.Polygon(pt_list)

def emf_tree(emf, tree,plane='xy',
              diameter_scale=1.0, linewidth=1.2,
              color=None, alpha=1.0,scalingFactor = 100, enhanceBifurcationPlotting = True, DiameterScaling = 2.0):

    plane0, plane1 = _plane2col(plane)
    ## the third element in the list is a boolean indicating if the section is a bifurcation

    section_segment_list = [(section, segment, len(section.children)>=2)
                            for section in iter_sections(tree)
                            for segment in iter_segments(section)]
    colors = [_get_color(color, section.type) for section, _,_ in section_segment_list]
    segs = [((seg[0][plane0], seg[0][plane1]),
                (seg[1][plane0], seg[1][plane1]))
            for _, seg,_ in section_segment_list]

    linewidth = _get_linewidth(
        tree,
        diameter_scale=diameter_scale,
        linewidth=linewidth, DiameterScaling = DiameterScaling)

    segs = (np.array(segs)*scalingFactor).astype('int32')
    segs[:,:,1] = -segs[:,:,1]       
    color = TREE_COLOR.get(tree.type, (0x00,0x00,0x00))
    for seg, color, width in zip(segs, colors, linewidth):
        emf.BeginPath()
        pen = emf.CreatePen(pyemf.PS_SOLID, int(width*scalingFactor),color)
        emf.SelectObject(pen)
        emf.MoveTo(seg[0][0],seg[0][1])
        emf.LineTo(seg[1][0],seg[1][1])
        emf.EndPath()
        emf.CloseFigure()
        emf.StrokeAndFillPath()
    if enhanceBifurcationPlotting:
        emf_bifurcation_plotting(emf, section_segment_list, scalingFactor, color)

def emf_bifurcation_plotting(emf, section_segment_list, scalingFactor, color):
    """Enhance the plotting of bifurcations by plotting a polygon around the bifurcation point"""
    for bif_point, _, isbif in section_segment_list:
        if isbif:
            parent_points = (np.array(bif_point.points[-1,:])*scalingFactor).astype('int32')
            parent_points[1] = -parent_points[1]
            for child in bif_point.children:
                child_points = (np.array(child.points[1,:])*scalingFactor).astype('int32')
                child_points[1] = -child_points[1]
                emf_polygon(emf, parent_points,child_points, color)


def emf_ScaleBar(emf, scalebar, scalingFactor = 100, Yoffset = 200):
    emf.SetTextAlign(pyemf.TA_BOTTOM|pyemf.TA_LEFT) 
    emf.SetTextColor((0x00,0x00,0x00))
    font = emf.CreateFont( 10*scalingFactor, 0, 0, 0, pyemf.FW_BOLD, 0, 0, 0,
                        pyemf.ANSI_CHARSET, pyemf.OUT_DEFAULT_PRECIS,
                        pyemf.CLIP_DEFAULT_PRECIS, pyemf.DEFAULT_QUALITY,
                        pyemf.DEFAULT_PITCH | pyemf.FF_DONTCARE, "Arial")
    emf.BeginPath()
    Yoffset = Yoffset*scalingFactor
    pen = emf.CreatePen(pyemf.PS_SOLID, int(1*scalingFactor),(0x00,0x00,0x00), [pyemf.PS_ENDCAP_ROUND,pyemf.PS_JOIN_BEVEL])
    emf.SelectObject(pen)
    emf.MoveTo(0,Yoffset)
    emf.LineTo(int(scalebar*scalingFactor),Yoffset)
    emf.EndPath()
    emf.CloseFigure()
    emf.StrokeAndFillPath()
    emf.SelectObject( font )
    emf.TextOut(int(scalebar*scalingFactor/2),Yoffset,str(scalebar)+' um')

def emf_morph(emf, morph,
               neurite_type=NeuriteType.all,
               plane='xy', scalingFactor=100, contourOn=True,rotationContour=0, scalebar=100,
               DiameterScaling = 2.0):
    """Plots a 2D figure of the morphology, that contains a soma and the neurites.
        scalebar is in um
    """   
    if scalebar is not None:
        emf_ScaleBar(emf, scalebar, scalingFactor)
    if contourOn:
        emf_contour(emf, morph, scalingFactor,rotationContour=rotationContour)
    emf_soma(emf, morph.soma, plane=plane,  scalingFactor=scalingFactor)
    for neurite in iter_neurites(morph, filt=tree_type_checker(neurite_type)):
        emf_tree(emf, neurite,  plane=plane, scalingFactor=scalingFactor, DiameterScaling = DiameterScaling)

def emf_contour(emf, morph,scalingFactor=100,rotationContour=0, contour_linewidth=100):
    contour_color = TREE_COLOR.get(NeuriteType.undefined, (0x00,0x00,0x00))
    pen = emf.CreatePen(pyemf.PS_DASH, contour_linewidth,contour_color)
    emf.SelectObject(pen)
    for idx, x in enumerate(morph.markers):
        ps = np.array(x.points)[:,:2]
        if rotationContour != 0:
            ps = rotate_contour(ps[:,:2], rotationContour)
        ps = (ps*scalingFactor).astype('int32')
        ps[:,1] = -ps[:,1]
        ps = ps.tolist()
        emf.Polyline(ps)

def orthognal_points(p1,p2,x):
    """Returns the coordinates of two points orthogonal to the line formed by points p1 and p2 with distance x"""
    v = p2-p1
    v = v/np.sqrt(np.sum(v**2))
    v = np.array([-v[1],v[0]])
    p3 = p1 + x*v
    p4 = p1 - x*v
    return p3.astype('int32'),p4.astype('int32')

def emf_polygon(emf, p1,p2,color, scalingFactor = 100):
    """Plots a polygon"""
    if np.abs(p1[COLS.R]-p2[COLS.R])/scalingFactor < 0.05:
        return
    p1a, p1b = orthognal_points(p1[:2],p2[:2], p1[COLS.R]//2)
    p1c, p1d = orthognal_points(p2[:2],p1[:2], p2[COLS.R]//2)
    shiftVector = (p2-p1)[:2]
    shiftVector = p1[COLS.R]*0.85*shiftVector/np.sqrt(np.sum(shiftVector**2))
    p1a = p1a + shiftVector.astype('int32')
    p1b = p1b + shiftVector.astype('int32')
    pen = emf.CreatePen(pyemf.PS_SOLID, 1,color)
    emf.SelectObject(pen)
    brush=emf.CreateSolidBrush(color)
    emf.SelectObject(brush)
    emf.Polygon([p1a,p1b,p1c,p1d])

def rotate_contour(contour, angle):
    """Rotate contour by angle in degrees."""
    angle = np.deg2rad(-angle)
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])
    return np.dot(contour, rotation_matrix)

def draw_morphorlogy_emf(morph,saveName, contourOn=True,rotationContour=0, width=210, height=297, dpi=300,unit='mm',
scalebar=100, DiameterScaling = 2.0, isSave = True, neuriteColors = None): 
    '''
    scalebar is in um
    isRadius = True: the COLS.R is the radius of the neurite; False: the COLS.R is the diameter of the neurite
    '''
    if neuriteColors is not None: # if neutriteColors is not None, update TREE_COLOR
        for idx, k in enumerate(neuriteColors): 
            if idx >= 4:
                break
            TREE_COLOR[k] = tuple([int(c*255) for c in neuriteColors[k]]) ## emf prefer intergers
    emf=pyemf.EMF(width,height,dpi,unit)
    emf.SetBkMode(pyemf.TRANSPARENT)
    emf_morph(emf, morph,plane='xy',scalingFactor=100, contourOn=contourOn, 
    scalebar=scalebar,rotationContour=rotationContour,DiameterScaling = DiameterScaling)
    if isSave:
        ret=emf.save(saveName+".emf")
        print(f"{saveName} save returns {ret}")

if __name__ == '__main__':
    cellName = '66'
    fname = cellName+"_mod.ASC"
    saveName = cellName+"_mod"
    neuron = morphor_nm.load_morphology(fname, somaType = SomaType.SOMA_CYLINDERS)
    draw_morphorlogy_emf(neuron, saveName)