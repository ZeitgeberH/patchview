"""
Heka Patchmaster .dat file reader 
tree = Bundle(fileName) ## load a tree bundle
pulseTree = tree.pul
data = tree.data[group_id, series_id, sweep_ind, trace_ind] # Load data for this trace
    
# stimTree = tree.pgf
stim = tree.StimWave[group_id, series_id, sweep_ind]
"""
import numpy as np
import re, struct, collections
import os
import datetime


class Struct(object):
    """High-level wrapper around struct.Struct that makes it a bit easier to
    unpack large, nested structures.
    """

    field_info = None
    size_check = None
    _fields_parsed = None

    def __init__(self, data, endian="<"):
        """Read the structure from *data* and return an ordered dictionary of
        fields.

        *data* may be a string or file.
        *endian* may be '<' or '>'
        """
        field_info = self._field_info()
        if not isinstance(data, (str, bytes)):
            data = data.read(self._le_struct.size)
        if endian == "<":
            items = self._le_struct.unpack(data)
        elif endian == ">":
            items = self._be_struct.unpack(data)
        else:
            raise ValueError("Invalid endian: %s" % endian)

        fields = collections.OrderedDict()

        i = 0
        for name, fmt, func in field_info:
            # pull item(s) out of the list based on format string
            if len(fmt) == 1 or fmt[-1] == "s":
                item = items[i]
                i += 1
            else:
                n = int(fmt[:-1])
                item = items[i : i + n]
                i += n

            # try unpacking sub-structure
            if isinstance(func, tuple):
                substr, func = func
                item = substr(item, endian)

            # None here means the field should be omitted
            if func is None:
                continue
            # handle custom massaging function
            if func is not True:
                item = func(item)
            fields[name] = item
            setattr(self, name, item)

        self.fields = fields

    @classmethod
    def _field_info(cls):
        if cls._fields_parsed is not None:
            return cls._fields_parsed

        fmt = ""
        fields = []
        for items in cls.field_info:
            if len(items) == 3:
                name, ifmt, func = items
            else:
                name, ifmt = items
                func = True

            if isinstance(ifmt, type) and issubclass(ifmt, Struct):
                func = (
                    ifmt,
                    func,
                )  # instructs to unpack with sub-struct before calling function
                ifmt = "%ds" % ifmt.size()
            elif len(ifmt) > 1 and re.match(r"\d*[xcbB?hHiIlLqQfdspP]", ifmt) is None:
                raise TypeError('Unsupported format string "%s"' % ifmt)

            fields.append((name, ifmt, func))
            fmt += ifmt
        cls._le_struct = struct.Struct("<" + fmt)
        cls._be_struct = struct.Struct(">" + fmt)
        cls._fields_parsed = fields
        if cls.size_check is not None:
            #            print(cls._le_struct.size, cls.size_check)
            assert cls._le_struct.size == cls.size_check
        return fields

    @classmethod
    def size(cls):
        cls._field_info()
        return cls._le_struct.size

    @classmethod
    def array(cls, x):
        """Return a new StructArray class of length *x* and using this struct
        as the array item type.
        """
        return type(
            cls.__name__ + "[%d]" % x,
            (StructArray,),
            {"item_struct": cls, "array_size": x},
        )

    def __repr__(self, indent=0):
        indent_str = "    " * indent
        r = indent_str + "%s(\n" % self.__class__.__name__
        if not hasattr(self, "fields"):
            r = r[:-1] + "<initializing>)"
            return r
        for k, v in self.fields.items():
            if isinstance(v, Struct):
                r += indent_str + "    %s = %s\n" % (
                    k,
                    v.__repr__(indent=indent + 1).lstrip(),
                )
            else:
                r += indent_str + "    %s = %r\n" % (k, v)
        r += indent_str + ")"
        return r

    def get_fields(self):
        """Recursively convert struct fields+values to nested dictionaries."""
        fields = self.fields.copy()
        for k, v in fields.items():
            if isinstance(v, StructArray):
                fields[k] = [x.get_fields() for x in v.array]
            elif isinstance(v, Struct):
                fields[k] = v.get_fields()
        return fields


class StructArray(Struct):
    item_struct = None
    array_size = None

    def __init__(self, data, endian="<"):
        if not isinstance(data, (str, bytes)):
            data = data.read(self.size())
        items = []
        isize = self.item_struct.size()
        for i in range(self.array_size):
            d = data[:isize]
            data = data[isize:]
            items.append(self.item_struct(d, endian))
        self.array = items

    def __getitem__(self, i):
        return self.array[i]

    @classmethod
    def size(self):
        return self.item_struct.size() * self.array_size

    def __repr__(self, indent=0):
        r = "    " * indent + "%s(\n" % self.__class__.__name__
        for item in self.array:
            r += item.__repr__(indent=indent + 1) + ",\n"
        r += "    " * indent + ")"
        return r


def cstr(byt):
    """Convert C string bytes to python string."""
    try:
        ind = byt.index(b"\0")
    except ValueError:
        return byt
    return byt[:ind].decode("utf-8", errors="ignore")


class BundleItem(Struct):
    field_info = [
        ("Start", "i"),
        ("Length", "i"),
        ("Extension", "8s", cstr),
    ]
    size_check = 16


class BundleHeader(Struct):
    ## DAT2 version file
    field_info = [
        ("Signature", "8s", cstr),
        ("Version", "32s", cstr),
        ("Time", "d"),
        ("Items", "i"),
        ("IsLittleEndian", "12s"),
        ("BundleItems", BundleItem.array(12)),
    ]
    size_check = 256


class TreeNode(Struct):
    """Struct that also represents a node in a Pulse file tree."""

    def __init__(self, fh, pul, level=0):
        self.level = level
        self.children = []
        endian = pul.endian

        # The record structure in the file may differ from our expected structure
        # due to version differences, so we read the required number of bytes, and
        # then pad or truncate before unpacking the record. This will probably
        # result in corrupt data in some situations..
        realsize = pul.level_sizes[level]
        data = fh.read(realsize)
        structsize = self.size()

        diff = structsize - realsize
        if diff != 0:
            #            print("TreeNode realsize: %g, struct size: %g, size diff: %g" % (realsize, structsize, diff))
            if diff > 0:
                data = data + b"\0" * diff
            else:
                data = data[:structsize]

        # initialize struct data
        Struct.__init__(self, data, endian)

        # Next read the number of children
        nchild = struct.unpack(endian + "i", fh.read(4))[0]

        level += 1
        if level >= len(pul.rectypes):
            return
        child_rectype = pul.rectypes[level]
        for i in range(nchild):
            self.children.append(child_rectype(fh, pul, level))

    def __getitem__(self, i):
        return self.children[i]

    def __len__(self):
        return len(self.children)

    def __iter__(self):
        return self.children.__iter__()

    def __repr__(self, indent=0):
        # Return a string describing this structure
        ind = "    " * indent
        srep = Struct.__repr__(self, indent)[:-1]  # exclude final parenthese
        srep += ind + "    children = %d,\n" % len(self)
        # srep += ind + 'children = [\n'
        # for ch in self:
        # srep += ch.__repr__(indent=indent+1) + ',\n'
        srep += ind + ")"
        return srep


class TraceRecord(TreeNode):
    """ """

    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("TraceCount", "i"),
        ("Data", "i"),
        ("DataPoints", "i"),
        ("InternalSolution", "i"),
        ("AverageCount", "i"),
        ("LeakCount", "i"),
        ("LeakTraces", "i"),
        ("DataKind", "h"),
        ("Filler1", "h", None),
        ("RecordingMode", "c"),
        ("AmplIndex", "c"),
        ("DataFormat", "c"),
        ("DataAbscissa", "c"),
        ("DataScaler", "d"),
        ("TimeOffset", "d"),
        ("ZeroData", "d"),
        ("YUnit", "8s", cstr),
        ("XInterval", "d"),
        ("XStart", "d"),
        ("XUnit", "8s", cstr),
        ("YRange", "d"),
        ("YOffset", "d"),
        ("Bandwidth", "d"),
        ("PipetteResistance", "d"),
        ("CellPotential", "d"),
        ("SealResistance", "d"),
        ("CSlow", "d"),
        ("GSeries", "d"),
        ("RsValue", "d"),
        ("GLeak", "d"),
        ("MConductance", "d"),
        ("LinkDAChannel", "i"),
        ("ValidYrange", "c"),
        ("AdcMode", "c"),
        ("AdcChannel", "h"),
        ("Ymin", "d"),
        ("Ymax", "d"),
        ("SourceChannel", "i"),
        ("ExternalSolution", "i"),
        ("CM", "d"),
        ("GM", "d"),
        ("Phase", "d"),
        ("DataCRC", "i"),
        ("CRC", "i"),
        ("GS", "d"),
        ("SelfChannel", "i"),
        ### add new fields compared to original heka_reader
        (
            "InterleaveSizeS",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				= 292; (* INT32 *)
        (
            "InterleaveSkip",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				= 296; (* INT32 *)
        ("ImageIndex", "i"),  # 		= fread(fh, 1, 'int32=>int32');%				= 300; (* INT32 *)
        (
            "Markers",
            "10d",
        ),  # 	= fread(fh, 10, 'double=>double');%				= 304; (* ARRAY[0..9] OF LONGREAL *)
        (
            "SECM_X",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double');%				= 384; (* LONGREAL *)
        (
            "SECM_Y",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double');%				= 392; (* LONGREAL *)
        (
            "SECM_Z",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double');%				= 400; (* LONGREAL *)
        (
            "TrHolding",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double');%				= 408; (* LONGREAL *)
        (
            "TcEnumerator",
            "i",
        ),  # 	= fread(fh, 1, 'int32=>int32');%				= 416; (* INT32 *)
        ("XTrace", "i"),  # 		= fread(fh, 1, 'int32=>int32');%				= 420; (* INT32 *)
        (
            "IntSolValue",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double');%				= 424; (* LONGREAL *)
        (
            "ExtSolValue",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double');%				= 432; (* LONGREAL *)
        (
            "IntSolName",
            "32s",
            cstr,
        ),  # 		= deblank(fread(fh, 32, 'uint8=>char')');%      = 440; (* String32Size *)
        (
            "ExtSolName",
            "32s",
            cstr,
        ),  # 		= deblank(fread(fh, 32, 'uint8=>char')');%      = 472; (* String32Size *)
        (
            "DataPedestal",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double');%				= 504; (* LONGREAL *)
    ]
    size_check = 512


"""
tr.TrInterleaveSize		= fread(fh, 1, 'int32=>int32');%				= 292; (* INT32 *)
tr.TrInterleaveSkip		= fread(fh, 1, 'int32=>int32');%				= 296; (* INT32 *)
tr.TrImageIndex			= fread(fh, 1, 'int32=>int32');%				= 300; (* INT32 *)
tr.TrMarkers			= fread(fh, 10, 'double=>double');%				= 304; (* ARRAY[0..9] OF LONGREAL *)
tr.TrSECM_X				= fread(fh, 1, 'double=>double');%				= 384; (* LONGREAL *)
tr.TrSECM_Y				= fread(fh, 1, 'double=>double');%				= 392; (* LONGREAL *)
tr.TrSECM_Z				= fread(fh, 1, 'double=>double');%				= 400; (* LONGREAL *)
tr.TrTrHolding			= fread(fh, 1, 'double=>double');%				= 408; (* LONGREAL *)
tr.TrTcEnumerator		= fread(fh, 1, 'int32=>int32');%				= 416; (* INT32 *)
tr.TrXTrace				= fread(fh, 1, 'int32=>int32');%				= 420; (* INT32 *)
tr.TrIntSolValue		= fread(fh, 1, 'double=>double');%				= 424; (* LONGREAL *)
tr.TrExtSolValue		= fread(fh, 1, 'double=>double');%				= 432; (* LONGREAL *)
tr.TrIntSolName			= deblank(fread(fh, 32, 'uint8=>char')');%      = 440; (* String32Size *)
tr.TrExtSolName			= deblank(fread(fh, 32, 'uint8=>char')');%      = 472; (* String32Size *)
tr.TrDataPedestal		= fread(fh, 1, 'double=>double');%				= 504; (* LONGREAL *)
"""


class SweepRecord(TreeNode):  ## Done.  validated with matlab version
    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("AuxDataFileOffset", "i"),
        ("StimCount", "i"),
        ("SweepCount", "i"),
        ("Time", "d"),
        ("Timer", "d"),
        ("SwUserParams", "2d"),
        ("PipPressure", "d"),
        ("RMSNoise", "d"),
        ("Temperature", "d"),
        ("OldIntSol", "i"),
        ("OldExtSol", "i"),
        ("DigitalIn", "h"),
        ("SweepKind", "h"),
        # ('DigitalOut','h'),
        ("Filler1", "i", None),
        ("Markers", "4d"),
        ("Filler2", "i", None),
        ("CRC", "i"),
        ("SwHolding", "16d"),
        ("UserParamEx", "8d"),
    ]
    size_check = 352


class V9_SweepRecord(TreeNode):
    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("AuxDataFileOffset", "i"),
        ("StimCount", "i"),
        ("SweepCount", "i"),
        ("Time", "d"),
        ("Timer", "d"),
        ("SwUserParams", "4d"),
        ("Temperature", "d"),
        ("OldIntSol", "i"),
        ("OldExtSol", "i"),
        ("DigitalIn", "h"),
        ("SweepKind", "h"),
        ("DigitalOut", "h"),
        ("Filler1", "i", None),
        ("Markers", "4d"),
        ("Filler2", "i", None),
        ("CRC", "i"),
        ("SwHolding", "16d"),
    ]
    ## according to Matlab Heka, but it could be 290 or 294. dependent on veriion of Heka
    ##TODO :  need clear this part!


#    size_check = 288


class UserParamDescrType(Struct):
    field_info = [
        ("Name", "32s", cstr),
        ("Unit", "8s", cstr),
    ]
    size_check = 40


class AmplifierState(Struct):
    field_info = [
        ("StateVersion", "8s", cstr),
        ("RealCurrentGain", "d"),
        ("RealF2Bandwidth", "d"),
        ("F2Frequency", "d"),
        ("RsValue", "d"),
        ("RsFraction", "d"),
        ("GLeak", "d"),
        ("CFastAmp1", "d"),
        ("CFastAmp2", "d"),
        ("CFastTau", "d"),
        ("CSlow", "d"),
        ("GSeries", "d"),
        ("StimDacScale", "d"),
        ("CCStimScale", "d"),
        ("VHold", "d"),
        ("LastVHold", "d"),
        ("VpOffset", "d"),
        ("VLiquidJunction", "d"),
        ("CCIHold", "d"),
        ("CSlowStimVolts", "d"),
        ("CCTrackVHold", "d"),
        ("TimeoutLength", "d"),
        ("SearchDelay", "d"),
        ("MConductance", "d"),
        ("MCapacitance", "d"),
        ("SerialNumber", "8s", cstr),
        ("E9Boards", "h"),
        ("CSlowCycles", "h"),
        ("IMonAdc", "h"),
        ("VMonAdc", "h"),
        ("MuxAdc", "h"),
        ("TstDac", "h"),
        ("StimDac", "h"),
        ("StimDacOffset", "h"),
        ("MaxDigitalBit", "h"),
        ("SpareInt1", "h", None),
        ("SpareInt2", "h", None),
        ("SpareInt3", "h", None),
        ("AmplKind", "c"),
        ("IsEpc9N", "c"),
        ("ADBoard", "c"),
        ("BoardVersion", "c"),
        ("ActiveE9Board", "c"),
        ("Mode", "c"),
        ("Range", "c"),
        ("F2Response", "c"),
        ("RsOn", "c"),
        ("CSlowRange", "c"),
        ("CCRange", "c"),
        ("CCGain", "c"),
        ("CSlowToTstDac", "c"),
        ("StimPath", "c"),
        ("CCTrackTau", "c"),
        ("WasClipping", "c"),
        ("RepetitiveCSlow", "c"),
        ("LastCSlowRange", "c"),
        ("Locked", "c"),
        ("CanCCFast", "c"),
        ("CanLowCCRange", "c"),
        ("CanHighCCRange", "c"),
        ("CanCCTracking", "c"),
        ("HasVmonPath", "c"),
        ("HasNewCCMode", "c"),
        ("Selector", "c"),
        ("HoldInverted", "c"),
        ("AutoCFast", "c"),
        ("AutoCSlow", "c"),
        ("HasVmonX100", "c"),
        ("TestDacOn", "c"),
        ("QMuxAdcOn", "c"),
        ("RealImon1Bandwidth", "d"),
        ("StimScale", "d"),
        ("Gain", "c"),
        ("Filter1", "c"),
        ("StimFilterOn", "c"),
        ("RsSlow", "c"),
        ("Old1", "c"),
        ("CCCFastOn", "c"),
        ("CCFastSpeed", "c"),
        ("F2Source", "c"),
        ("TestRange", "c"),
        ("TestDacPath", "c"),
        ("MuxChannel", "c"),
        ("MuxGain64", "c"),
        ("VmonX100", "c"),
        ("IsQuadro", "c"),
        ("SpareBool4", "c", None),
        ("SpareBool5", "c", None),
        ("StimFilterHz", "d"),
        ("RsTau", "d"),
        ("FilterOffsetDac", "h"),
        ("ReferenceDac", "h"),
        ("SpareInt6", "h", None),
        ("SpareInt7", "h", None),
        ("Spares1", "24s", None),
        ("CalibDate", "16s"),
        ("SelHold", "d"),
        ("Spares2", "32s", None),
    ]
    size_check = 400


class LockInParams(Struct):
    field_info = [
        ("ExtCalPhase", "d"),
        ("ExtCalAtten", "d"),
        ("PLPhase", "d"),
        ("PLPhaseY1", "d"),
        ("PLPhaseY2", "d"),
        ("UsedPhaseShift", "d"),
        ("UsedAttenuation", "d"),
        ("Spares2", "8s", None),
        ("ExtCalValid", "?"),
        ("PLPhaseValid", "?"),
        ("LockInMode", "c"),
        ("CalMode", "c"),
        ("Spares", "28s", None),
    ]
    size_check = 96


class V9_SeriesRecord(TreeNode):  ## Done!
    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("Comment", "80s", cstr),
        ("SeriesCount", "i"),
        ("NumberSweeps", "i"),
        ("AmplStateOffset", "i"),
        ("AmplStateSeries", "i"),
        ("MethodTag", "i"),
        ("Time", "d"),
        ("PageWidth", "d"),
        ("SwUserParamDescr", UserParamDescrType.array(4)),
        ("MethodName", "32s", cstr),
        ("SeUserParams", "4d"),
        ("LockInParams", LockInParams),
        ("AmplifierState", AmplifierState),
        ("Username", "80s", cstr),
        ("UserParamDescr", UserParamDescrType.array(4)),
        ("Filler1", "i", None),
        ("CRC", "i"),
        ("UserParams2", "4d"),
        ("UserParamDescr2", UserParamDescrType.array(4)),
        ("ScanParams", "96c"),  ## 96 uint8
    ]
    size_check = 1408


class SeriesRecord(TreeNode):  ## Done! validated with Matlab HEKA importer!
    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("Comment", "80s", cstr),
        ("SeriesCount", "i"),
        ("NumberSweeps", "i"),
        ("AmplStateFlag", "i"),
        ("AmplStateRef", "i"),
        ("MethodTag", "i"),
        ("Time", "d"),
        ("PageWidth", "d"),
        ("SwUserParamDescr", UserParamDescrType.array(2)),
        ("Filler1", "80s", None),
        ("MethodName", "32s", cstr),
        ("PhotoParams1", "4d"),
        ("OldLockInParams", LockInParams),
        ("OldAmpState", AmplifierState),
        ("Username", "80s", cstr),
        ("PhotoParams2", UserParamDescrType.array(4)),
        ("Filler1", "i", None),
        ("CRC", "i"),
        ("UserParams2", "4d"),
        ("UserParamDescr2", UserParamDescrType.array(4)),
        ("ScanParams", "96c"),  ## 96 uint8
        ("UserDescr2", UserParamDescrType.array(8)),
    ]
    size_check = 1728


class GroupRecord(TreeNode):

    field_info = [
        ("Mark", "i"),
        ("Label", "32s", cstr),
        ("Text", "80s", cstr),
        ("ExperimentNumber", "i"),
        ("GroupCount", "i"),
        ("CRC", "i"),
        ("MatrixWidth", "d"),
        ("MatrixHeight", "d"),
    ]
    size_check = 144


class Pulsed(TreeNode):
    field_info = [
        ("Version", "i"),
        ("Mark", "i"),
        ("VersionName", "32s", cstr),
        ("AuxFileName", "80s", cstr),
        ("RootText", "400s", cstr),
        ("StartTime", "d"),
        ("MaxSamples", "i"),
        ("CRC", "i"),
        ("Features", "h"),
        ("Filler1", "h", None),
        ("Filler2", "i", None),
        ("TcEnumerator", "32h"),
        ("TcKind", "32c"),
    ]
    size_check = 640

    rectypes = [None, GroupRecord, SeriesRecord, SweepRecord, TraceRecord]

    def __init__(self, bundle, offset=0, size=None):
        fh = open(bundle.file_name, "rb")
        # pdb.set_trace()
        fh.seek(offset)

        # read .pul header
        magic = fh.read(4)
        if magic == b"eerT":
            self.endian = "<"
        elif magic == b"Tree":
            self.endian = ">"
        elif magic == b"DAT1":
            self.endian = ">"
        else:
            raise RuntimeError("Bad file magic: %s" % magic)

        levels = struct.unpack(self.endian + "i", fh.read(4))[0]
        # read size of each level (one int per level)
        self.level_sizes = []
        for i in range(levels):
            size = struct.unpack(self.endian + "i", fh.read(4))[0]
            self.level_sizes.append(size)

        TreeNode.__init__(self, fh, self)


class Pulsed9(TreeNode):
    field_info = [
        ("Version", "i"),
        ("Mark", "i"),
        ("VersionName", "32s", cstr),
        ("AuxFileName", "80s", cstr),
        ("RootText", "400s", cstr),
        ("StartTime", "d"),
        ("MaxSamples", "i"),
        ("CRC", "i"),
        ("Features", "h"),
        ("Filler1", "h", None),
        ("Filler2", "i", None),
    ]
    size_check = 544

    rectypes = [
        None,
        GroupRecord,  ## no changes in group record between version 9 and version 1000
        V9_SeriesRecord,
        V9_SweepRecord,
        TraceRecord,  ## no changes in tracerecord between version 9 and version 1000
    ]

    def __init__(self, bundle, offset=0, size=None):
        fh = open(bundle.file_name, "rb")
        fh.seek(offset)

        # read .pul header
        magic = fh.read(4)
        if magic == b"eerT":
            self.endian = "<"
        elif magic == b"Tree":
            self.endian = ">"
        else:
            raise RuntimeError("Bad file magic: %s" % magic)

        levels = struct.unpack(self.endian + "i", fh.read(4))[0]

        # read size of each level (one int per level)
        self.level_sizes = []
        for i in range(levels):
            size = struct.unpack(self.endian + "i", fh.read(4))[0]
            self.level_sizes.append(size)

        TreeNode.__init__(self, fh, self)


class Data(object):
    def __init__(self, bundle, offset=0, size=None):
        self.bundle = bundle
        self.offset = offset

    def __getitem__(self, *args):
        index = args[0]
        assert len(index) == 4
        pul = self.bundle.pul
        trace = pul[index[0]][index[1]][index[2]][index[3]]
        fh = open(self.bundle.file_name, "rb")
        fh.seek(trace.Data)
        fmt = bytearray(trace.DataFormat)[0]
        dtype = [np.int16, np.int32, np.float32, np.float64][fmt]
        dByte = [2, 4, 4, 8][fmt]
        nItemsPerBlock = np.int32(trace.InterleaveSizeS / dByte)
        TotalBytes = trace.DataPoints * dByte
        # print('{:f}, {:f}'.format(trace.DataPoints, trace.InterleaveSizeS))
        if trace.DataPoints >= trace.InterleaveSizeS and trace.InterleaveSizeS != 0:
            print("long block")
            ### there is a mixture of data points (count) and bytes!
            data = np.fromfile(fh, count=nItemsPerBlock, dtype=dtype)
            dataRead = trace.InterleaveSizeS
            data2Read = TotalBytes - dataRead  ## in bytes
            c = 0
            while data2Read > 0:
                fh.seek(
                    trace.InterleaveSkip - trace.InterleaveSizeS, os.SEEK_CUR
                )  ## skip the skip-block
                c = c + 1
                #                print(c)
                if data2Read < trace.InterleaveSizeS:  ## less than a block size
                    data0 = np.fromfile(
                        fh, count=np.int(data2Read / dByte), dtype=dtype
                    )
                    data = np.concatenate((data, data0))
                    break
                else:  ## larger than a block size
                    data0 = np.fromfile(fh, count=nItemsPerBlock, dtype=dtype)
                    data = np.concatenate((data, data0))
                    dataRead = trace.InterleaveSizeS + dataRead
                    data2Read = TotalBytes - dataRead
        else:
            data = np.fromfile(fh, count=trace.DataPoints, dtype=dtype)
        return data * trace.DataScaler + trace.ZeroData


class StimulationRecord(TreeNode):
    """
    (* StimulationRecord = RECORD *)
    """

    ### Long real: d
    ### Byte:  c
    field_info = [
        ("Mark", "i"),  #               =   0; (* INT32 *)
        ("EntryName", "32s", cstr),  #          =   4; (* String32Type *)
        ("FileName", "32s", cstr),  #           =  36; (* String32Type *)
        ("AnalName", "32s", cstr),  #           =  68; (* String32Type *)
        ("DataStartSegment", "i"),  #    = 100; (* INT32 *)
        ("DataStartTime", "d"),  #      = 104; (* LONGREAL *)
        ("SampleInterval", "d"),  #     = 112; (* LONGREAL *)
        ("SweepInterval", "d"),  #      = 120; (* LONGREAL *)
        ("LeakDelay", "d"),  #         = 128; (* LONGREAL *)
        ("FilterFactor", "d"),  #       = 136; (* LONGREAL *)
        ("NumberSweeps", "i"),  #       = 144; (* INT32 *)
        ("NumberLeaks", "i"),  #        = 148; (* INT32 *)
        ("NumberAverages ", "i"),  #     = 152; (* INT32 *)
        ("ActualAdcChannels", "i"),  #   = 156; (* INT32 *)
        ("ActualDacChannels ", "i"),  #  = 160; (* INT32 *)
        ("ExtTrigger", "c"),  #         = 164; (* BYTE *)
        ("NoStartWait", "h"),  #        = 165; (* BOOLEAN *)
        ("UseScanRates", "h"),  #       = 166; (* BOOLEAN *)
        ("NoContAq", "h"),  #           = 167; (* BOOLEAN *)
        ("HasLockIn", "h"),  #          = 168; (* BOOLEAN *)
        ("OldStartMacKind", "c"),  # = 169; (* CHAR *)
        ("OldEndMacKind", "h"),  #   = 170; (* BOOLEAN *)
        ("AutoRange", "c"),  #          = 171; (* BYTE *)
        ("BreakNext", "h"),  #          = 172; (* BOOLEAN *)
        ("IsExpanded", "h"),  #         = 173; (* BOOLEAN *)
        ("LeakCompMode", "h"),  #       = 174; (* BOOLEAN *)
        ("HasChirp", "h"),  #            = 175; (* BOOLEAN *)
        ("OldStartMacro", "32s", cstr),  #   = 176; (* String32Type *)
        ("OldEndMacro", "32s", cstr),  #     = 208; (* String32Type *)
        ("IsGapFree", "h"),  #           = 240; (* BOOLEAN *)
        ("HandledExternally ", "h"),  #   = 241; (* BOOLEAN *)
        ("Filler1", "i"),  #         = 242; (* BOOLEAN *)
        ("Filler2", "i"),  #         = 243; (* BOOLEAN *)
        ("CRC", "i"),  #                = 244; (* CARD32 *)
    ]


#    size_check = 248


class ChannelRecord(TreeNode):
    """
    set fileds of Channel record
    """

    field_info = [
        ("Mark", "i"),
        ("LinkedChannel", "i"),
        (
            "CompressionFactor",
            "i",
        ),  # 	= fread(fh, 1, 'int32=>int32');%				=   8; (* INT32 *)
        (
            "YUnit",
            "8s",
            cstr,
        ),  # 			= deblank(fread(fh, 8, 'uint8=>char')');%       =  12; (* String8Type *)
        ("AdcChannel", "h"),  # 		= fread(fh, 1, 'int16=>int16');%				=  20; (* INT16 *)
        ("AdcMode", "c"),  # 			= fread(fh, 1, 'uint8=>uint8');%				=  22; (* BYTE *)
        (
            "DoWrite",
            "b",
        ),  # 			= fread(fh, 1, 'uint8=>logical');%           	=  23; (* BOOLEAN *)
        (
            "stLeakStore",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>uint8');%				=  24; (* BYTE *)
        ("AmplMode", "c"),  # 		= fread(fh, 1, 'uint8=>uint8');%				=  25; (* BYTE *)
        (
            "OwnSegTime",
            "c",
        ),  # 	= fread(fh, 1, 'uint8=>logical');%				=  26; (* BOOLEAN *)
        (
            "SetLastSegVmemb",
            "c",
        ),  # 	= fread(fh, 1, 'uint8=>logical');%				=  27; (* BOOLEAN *)
        ("DacChannel", "h"),  # 		= fread(fh, 1, 'int16=>int16');%				=  28; (* INT16 *)
        ("DacMode", "c"),  # 			= fread(fh, 1, 'uint8=>uint8');%				=  30; (* BYTE *)
        (
            "HasLockInSquare",
            "c",
        ),  # 		= fread(fh, 1, 'uint8=>uint8');%				=  31; (* BYTE *)
        (
            "RelevantXSegment",
            "i",
        ),  # 	= fread(fh, 1, 'int32=>int32');%				=  32; (* INT32 *)
        (
            "RelevantYSegment",
            "i",
        ),  # = fread(fh, 1, 'int32=>int32');%				=  36; (* INT32 *)
        (
            "DacUnit",
            "8s",
            cstr,
        ),  # 		= deblank(fread(fh, 8, 'uint8=>char')');%       =  40; (* String8Type *)
        (
            "Holding",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				=  48; (* LONGREAL *)
        (
            "LeakHolding",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				=  56; (* LONGREAL *)
        (
            "LeakSize",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double') ;%				=  64; (* LONGREAL *)
        ("LeakHoldMode", "c"),  # 	= fread(fh, 1, 'uint8=>uint8');%				=  72; (* BYTE *)
        (
            "LeakAlternate",
            "c",
        ),  # 	= fread(fh, 1, 'uint8=>logical');%				=  73; (* BOOLEAN *)
        (
            "AltLeakAveraging",
            "c",
        ),  # = fread(fh, 1, 'uint8=>logical');%				=  74; (* BOOLEAN *)
        (
            "LeakPulseOn",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>logical');%				=  75; (* BOOLEAN *)
        (
            "StimToDacID",
            "h",
        ),  # 			= fread(fh, 1, 'int16=>int16');%				=  76; (* SET16 *)
        (
            "CompressionMode",
            "h",
        ),  # 		= fread(fh, 1, 'int16=>int16');%				=  78; (* SET16 *)
        (
            "CompressionSkip",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				=  80; (* INT32 *)
        (
            "DacBit",
            "h",
        ),  # 				= fread(fh, 1, 'int16=>int16');%            	=  84; (* INT16 *)
        (
            "HasLockInSine",
            "c",
        ),  # 		= fread(fh, 1, 'uint8=>logical');%				=  86; (* BOOLEAN *)
        ("BreakMode", "c"),  # 			= fread(fh, 1, 'uint8=>uint8');%				=  87; (* BYTE *)
        ("ZeroSeg", "i"),  # 				= fread(fh, 1, 'int32=>int32');%				=  88; (* INT32 *)
        ("StimSweep", "i"),  # 			= fread(fh, 1, 'int32=>int32');%				=  92; (* INT32 *)
        (
            "Sine_Cycle",
            "d",
        ),  # 			= fread(fh, 1, 'double=>double') ;%				=  96; (* LONGREAL *)
        (
            "Sine_Amplitude",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 104; (* LONGREAL *)
        (
            "LockIn_VReversal",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double') ;%				= 112; (* LONGREAL *)
        (
            "Chirp_StartFreq",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 120; (* LONGREAL *)
        (
            "Chirp_EndFreq",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 128; (* LONGREAL *)
        (
            "Chirp_MinPoints",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 136; (* LONGREAL *)
        (
            "Square_NegAmpl",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 144; (* LONGREAL *)
        (
            "Square_DurFactor",
            "d",
        ),  # 	= fread(fh, 1, 'double=>double') ;%				= 152; (* LONGREAL *)
        (
            "LockIn_Skip",
            "i",
        ),  # 			= fread(fh, 1, 'int32=>int32');%				= 160; (* INT32 *)
        (
            "Photo_MaxCycles",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				= 164; (* INT32 *)
        (
            "Photo_SegmentNo",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				= 168; (* INT32 *)
        (
            "LockIn_AvgCycles",
            "i",
        ),  # 	= fread(fh, 1, 'int32=>int32');%				= 172; (* INT32 *)
        (
            "Imaging_RoiNo",
            "i",
        ),  # 		= fread(fh, 1, 'int32=>int32');%				= 176; (* INT32 *)
        (
            "Chirp_Skip",
            "i",
        ),  # 			= fread(fh, 1, 'int32=>int32');%				= 180; (* INT32 *)
        (
            "Chirp_Amplitude",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double');%				= 184; (* LONGREAL *)
        (
            "Photo_Adapt",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>uint8');%				= 192; (* BYTE *)
        ("Sine_Kind", "c"),  # 			= fread(fh, 1, 'uint8=>uint8');%				= 193; (* BYTE *)
        (
            "Chirp_PreChirp",
            "c",
        ),  # 		= fread(fh, 1, 'uint8=>uint8');%				= 194; (* BYTE *)
        (
            "Sine_Source",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>uint8');%				= 195; (* BYTE *)
        (
            "Square_NegSource",
            "c",
        ),  # 	= fread(fh, 1, 'uint8=>uint8');%				= 196; (* BYTE *)
        (
            "Square_PosSource",
            "c",
        ),  # 	= fread(fh, 1, 'uint8=>uint8');%				= 197; (* BYTE *)
        ("Chirp_Kind", "c"),  # 			= fread(fh, 1, 'uint8=>uint8');%				= 198; (* BYTE *)
        (
            "Chirp_Source",
            "c",
        ),  # 		= fread(fh, 1, 'uint8=>uint8');%				= 199; (* BYTE *)
        (
            "DacOffset",
            "d",
        ),  # 			= fread(fh, 1, 'double=>double') ;%				= 200; (* LONGREAL *)
        (
            "AdcOffset",
            "d",
        ),  # 			= fread(fh, 1, 'double=>double') ;%				= 208; (* LONGREAL *)
        (
            "TraceMathFormat",
            "c",
        ),  # 		= fread(fh, 1, 'uint8=>uint8');%				= 216; (* BYTE *)
        (
            "HasChirp",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>logical');%				= 217; (* BOOLEAN *)
        (
            "Square_Kind",
            "c",
        ),  # 			= fread(fh, 1, 'uint8=>uint8');%				= 218; (* BYTE *)
        (
            "Filler1",
            "5s",
            cstr,
        ),  # 			= fread(fh, 5, 'uint8=>char');%					= 219; (* ARRAY[0..5] OF CHAR *)
        (
            "Square_BaseIncr",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double');%				= 224; (* LONGREAL *)
        (
            "Square_Cycle",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 232; (* LONGREAL *)
        (
            "Square_PosAmpl",
            "d",
        ),  # 		= fread(fh, 1, 'double=>double') ;%				= 240; (* LONGREAL *)
        (
            "CompressionOffset",
            "i",
        ),  # 	= fread(fh, 1, 'int32=>int32');%				= 248; (* INT32 *)
        ("PhotoMode", "i"),  # 			= fread(fh, 1, 'int32=>int32');%				= 252; (* INT32 *)
        (
            "BreakLevel",
            "d",
        ),  # 			= fread(fh, 1, 'double=>double') ;%				= 256; (* LONGREAL *)
        (
            "TraceMath",
            "128s",
            cstr,
        ),  # 			= deblank(fread(fh,128,'uint8=>char')');%       = 264; (* String128Type *)
        #            ('OldCRC','i'),		#				= fread(fh, 1, 'int32=>int32');%                = 268; (* CARD32 *)
        ("Filler2", "i"),  # 			= fread(fh, 1, 'int32=>int32');%				= 392; (* INT32 *)
        (
            "CRC",
            "i",
        ),  # 					= fread(fh, 1, 'int32=>int32');%                = 396; (* CARD32 *)
    ]
    size_check = 400


#%%
class StimSegmentRecord(TreeNode):
    """
    (* StimSegmentRecord = RECORD *)
    """

    field_info = [
        ("Mark", "i"),  #               =   0; (* INT32 *)
        ("Class", "c"),  #              =   4; (* BYTE *)
        ("StoreKind", "c"),  #      =   5; (* BYTE *)
        ("VoltageIncMode", "c"),  #     =   6; (* BYTE *)
        ("DurationIncMode", "c"),  #  =   7; (* BYTE *)
        ("Voltage", "d"),  #    =   8; (* LONGREAL *)
        ("VoltageSource", "i"),  #      =  16; (* INT32 *)
        ("DeltaVFactor", "d"),  #       =  20; (* LONGREAL *)
        ("DeltaVIncrement", "d"),  #   =  28; (* LONGREAL *)
        ("Duration", "d"),  #      =  36; (* LONGREAL *)
        ("DurationSource", "i"),  #    =  44; (* INT32 *)
        ("DeltaTFactor", "d"),  #     =  48; (* LONGREAL *)
        ("DeltaTIncrement", "d"),  #   =  56; (* LONGREAL *)
        ("Filler1", "i"),  #        =  64; (* INT32 *)
        ("CRC", "i"),  #             =  68; (* CARD32 *)
        ("ScanRate", "d"),  #        =  72; (* LONGREAL *)
    ]
    size_check = 80


#%%
class StimTree(TreeNode):
    field_info = [
        ("Version", "i"),
        ("Mark", "i"),
        ("VersionName", "32s", cstr),
        ("MaxSamples", "i"),
        ("Filler1", "i"),
        ("Params", "10d"),
        ("ParamText", "320s", cstr),
        ("Reserved", "32i"),
        ("Filler2", "i"),
        ("CRC", "i"),
    ]
    size_check = 584

    rectypes = [
        None,
        StimulationRecord,
        ChannelRecord,  ## not used for now. ONly the second level has useful information. extract data from this level!
        StimSegmentRecord,
    ]

    def __init__(self, bundle, offset=0, size=None):
        fh = open(bundle.file_name, "rb")
        fh.seek(offset)

        # read .stim header
        stimHeaderIdx = 4
        magic = fh.read(stimHeaderIdx)
        if magic == b"eerT":
            self.endian = "<"
        elif magic == b"Tree":
            self.endian = ">"
        else:
            raise RuntimeError("Bad file magic: %s" % magic)

        levels = struct.unpack(self.endian + "i", fh.read(stimHeaderIdx))[0]

        # read size of each level (one int per level)
        self.level_sizes = []
        for i in range(levels):
            size = struct.unpack(self.endian + "i", fh.read(stimHeaderIdx))[0]
            self.level_sizes.append(size)

        TreeNode.__init__(self, fh, self)


def datTime2RealTime(tc):
    """Convert from timestamp to human-readable time"""
    tc = tc - 1580970496
    if tc < 0:
        tc = tc + 4294967296
    tc = tc + 9561652096
    T_hsa = datetime.datetime.utcfromtimestamp(tc) - datetime.timedelta(134774)
    return T_hsa


class Bundle(object):

    item_classes = {
        ".pul": Pulsed,  ## pulse tree for file version 1000
        ".pul9": Pulsed9,  ## pulse tree version 9
        ".dat": Data,
        ".pgf": StimTree,
    }

    def __init__(self, file_name):
        self.file_name = file_name
        fh = open(file_name, "rb")
        self.filehandle = fh
        # Read header assuming little endiam
        endian = "<"
        self.header = BundleHeader(fh, endian)
        #        print(self.header)
        # print("Header version: " + self.header.Version)
        self.session_start_time = datTime2RealTime(self.header.Time)
        self.signature = self.header.Signature
        if self.signature == "DAT1":
            self.isBundled = False
            print("WARNING: this dat file is not bundled!")
            print("Need to read .pul file for data tree!")
        else:
            self.isBundled = True
        # print("Recording start time:"+ str(self.session_start_time))
        if (
            self.header.Version[:7] == "v2x90.2"
        ):  ## can not handle files older than v2x90.2
            self.fileVersion = 9
        else:
            self.fileVersion = 1000

        self.endian = "little"
        # print('header size', self.header.size())
        # If the header is bad, re-read using big endian
        if self.header.IsLittleEndian[0] == b"\0":
            endian = ">"
            fh.seek(0)
            self.header = BundleHeader(fh, endian)
            self.endian = "big"
        self.catalog = {}
        for item in self.header.BundleItems:
            item.instance = None
            ext = item.Extension
            self.catalog[ext] = item

    def closeDatFile(self):
        """
        close opened file
        """
        return self.filehandle.close()

    @property
    def pul(self):
        """The Pulsed object from this pul."""
        if self.isBundled:
            return self._get_item_instance(".pul")
        else:
            return Pulsed(self, offset=0)

    @property
    def data(self):
        """The Data object from this bundle."""
        return self._get_item_instance(".dat")

    @property
    def pgf(self):
        """The pgf object from this bundle."""
        return self._get_item_instance(".pgf")

    def _get_item_instance(self, ext):
        # pdb.set_trace()
        if ext not in self.catalog:
            return None
        item = self.catalog[ext]
        if item.instance is None:
            if ext == ".pul" and self.fileVersion == 9:
                cls = self.item_classes[".pul9"]
            elif ext == ".pul" and self.fileVersion == 10:
                cls = self.item_classes[".pul"]
            else:
                cls = self.item_classes[ext]
            item.instance = cls(self, item.Start, item.Length)
        return item.instance

    def stim(self, index):
        """
        Pgf@stimulation == Pul@series. All stimulation are collapase into the pgf tree,
        so if there are multiple groups in Pulse tree, then the series index
        need to be accumulated from the first group, rather than from that group
        level.
        """
        stimTree = self.pgf
        assert len(index) >= 2  ## Do this at channel level
        groupID = index[0]  ## Group index.
        stimID = 0  ## start at the first stimulaition record
        if (
            groupID >= 1
        ):  ## if this is not the first group, need to check accumated series index
            for g in range(groupID - 1):
                stimID = stimID + len(
                    self.pul[g].children
                )  ## get number of series in this group
        stimID = stimID + index[1]  ## plus current series index
        stimulationRecord = stimTree.children[stimID]
        stimulatingElectrode = stimulationRecord.EntryName
        ## index for channel record. We just use the first channel.
        # Rest channes within same simulation record is the same
        iC = 0
        chanRecord = stimulationRecord.children[iC]
        #        DacMode = int.from_bytes(chanRecord.DacMode, byteorder =  self.endian)  ## relative or abolute
        Stim2DacID = chanRecord.StimToDacID
        #        stim2dac = ['StimScale','StimScale, Relative','FileTemplate','FileTemplate, Relative']
        if Stim2DacID == 1:
            stim2dacType = "StimScale"
        elif Stim2DacID == 3:
            stim2dacType = "StimScale, Relative"
        else:
            print("unknow stimToDacID. Assume relative scale")
            Stim2DacID = 1
        Holding = chanRecord.Holding * 1000.0  ## pA
        sampleInteval = stimulationRecord.SampleInterval
        stim = []
        cumSamples = 0  ## cumulated stimulation time
        stimInfo = []
        if len(index) == 4:
            traceRecord = self.pul[index[0]][index[1]][index[2]][index[3]]
            #            print(traceRecord)
            sealResistance = traceRecord.SealResistance
        else:
            sealResistance = np.nan
        iSweep = index[2]  ## index for sweep. This is only used to generate the
        nsweep = len(chanRecord.children)
        ## the current value
        for idx, seg in enumerate(
            chanRecord.children
        ):  ## go through segment one by one
            segSamples = int(seg.Duration / sampleInteval)
            #            segV0 = seg.Voltage  ## initial voltage
            # if Stim2DacID == 3:  ## relative mode
            #     seg_V = 1000.0*(seg.Voltage + iSweep*seg.DeltaVIncrement*seg.DeltaVFactor) +Holding
            # else:  ## abosolute
            if idx == 0 or idx == nsweep - 1:
                seg_V = Holding
            else:
                if Stim2DacID != 3:
                    seg_V = 1000.0 * (
                        seg.Voltage + iSweep * seg.DeltaVIncrement * seg.DeltaVFactor
                    )
                else:
                    seg_V = (
                        1000.0
                        * (
                            seg.Voltage
                            + iSweep * seg.DeltaVIncrement * seg.DeltaVFactor
                        )
                        + Holding
                    )

            #            seg_dT = seg.DeltaTIncrement*seg.DeltaTFactor ## probably not necessary for stepwise type of stimuli
            stim.extend(np.ones((segSamples,)) * seg_V)
            stim_ = {
                "EntryName": stimulatingElectrode,
                "start": cumSamples,
                "end": cumSamples + segSamples,
                "duration": seg.Duration,
                "amplitude": seg_V,
                "sampleInteval": sampleInteval,
                "Vholding": Holding,
                "StimToDac": stim2dacType,
                "sealresistance": sealResistance,
            }
            cumSamples = cumSamples + segSamples
            stimInfo.append(stim_)

        time = np.arange(cumSamples) * sampleInteval
        #        stimInfo = pd.DataFrame(stimInfo)
        return time, stim, stimInfo

    def __repr__(self):
        return "Bundle(%r)" % list(self.catalog.keys())
