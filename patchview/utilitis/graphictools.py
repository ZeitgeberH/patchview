def diameterToMPLPixels(dpi, scalingFactor=1000):
    """scaling factor for a physical diameter in um to the unit of matplotlib pixels.
    Matplotlib figures use Points per inch (ppi) of 72, so 1 inch = 72 points.
    In Illustrator, high resolution is 300 ppi, low resolution is 72 ppi.
    A4 size paper has 841.89 pts (210 mm)  width and 595.28 pts (297 mm) height
    1 point = 1/72 inch = 0.352777778 mm = 0.0138888889 inch
    Most elements like lines, markers, texts have a size given in points.
    The default DPI in 100
    https://stackoverflow.com/questions/47633546/relationship-between-dpi-and-figure-size
    scalingFactor is to magnify the size
    1 mm = 0.0393701 inch
    """
    return scalingFactor*0.0393701*dpi/72