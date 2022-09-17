"""Legacy module, replaced by neurom.features"""
from warnings import warn
from patchview.neurom.features import neuritefunc as _neuritefunc
from patchview.neurom.features import neuronfunc as _neuronfunc
from patchview.neurom.features import sectionfunc

from patchview.neurom.features import (
    NEURITEFEATURES,
    NEURONFEATURES,
    get,
    register_neurite_feature,
)
