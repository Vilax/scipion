# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.35
#
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _XmippData
import new
new_instancemethod = new.instancemethod
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class PySwigIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PySwigIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PySwigIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError, "No constructor defined"
    __repr__ = _swig_repr
    __swig_destroy__ = _XmippData.delete_PySwigIterator
    __del__ = lambda self : None;
    def value(*args): return _XmippData.PySwigIterator_value(*args)
    def incr(*args): return _XmippData.PySwigIterator_incr(*args)
    def decr(*args): return _XmippData.PySwigIterator_decr(*args)
    def distance(*args): return _XmippData.PySwigIterator_distance(*args)
    def equal(*args): return _XmippData.PySwigIterator_equal(*args)
    def copy(*args): return _XmippData.PySwigIterator_copy(*args)
    def next(*args): return _XmippData.PySwigIterator_next(*args)
    def previous(*args): return _XmippData.PySwigIterator_previous(*args)
    def advance(*args): return _XmippData.PySwigIterator_advance(*args)
    def __eq__(*args): return _XmippData.PySwigIterator___eq__(*args)
    def __ne__(*args): return _XmippData.PySwigIterator___ne__(*args)
    def __iadd__(*args): return _XmippData.PySwigIterator___iadd__(*args)
    def __isub__(*args): return _XmippData.PySwigIterator___isub__(*args)
    def __add__(*args): return _XmippData.PySwigIterator___add__(*args)
    def __sub__(*args): return _XmippData.PySwigIterator___sub__(*args)
    def __iter__(self): return self
PySwigIterator_swigregister = _XmippData.PySwigIterator_swigregister
PySwigIterator_swigregister(PySwigIterator)

class vectori(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectori, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectori, name)
    __repr__ = _swig_repr
    def iterator(*args): return _XmippData.vectori_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _XmippData.vectori___nonzero__(*args)
    def __len__(*args): return _XmippData.vectori___len__(*args)
    def pop(*args): return _XmippData.vectori_pop(*args)
    def __getslice__(*args): return _XmippData.vectori___getslice__(*args)
    def __setslice__(*args): return _XmippData.vectori___setslice__(*args)
    def __delslice__(*args): return _XmippData.vectori___delslice__(*args)
    def __delitem__(*args): return _XmippData.vectori___delitem__(*args)
    def __getitem__(*args): return _XmippData.vectori___getitem__(*args)
    def __setitem__(*args): return _XmippData.vectori___setitem__(*args)
    def append(*args): return _XmippData.vectori_append(*args)
    def empty(*args): return _XmippData.vectori_empty(*args)
    def size(*args): return _XmippData.vectori_size(*args)
    def clear(*args): return _XmippData.vectori_clear(*args)
    def swap(*args): return _XmippData.vectori_swap(*args)
    def get_allocator(*args): return _XmippData.vectori_get_allocator(*args)
    def begin(*args): return _XmippData.vectori_begin(*args)
    def end(*args): return _XmippData.vectori_end(*args)
    def rbegin(*args): return _XmippData.vectori_rbegin(*args)
    def rend(*args): return _XmippData.vectori_rend(*args)
    def pop_back(*args): return _XmippData.vectori_pop_back(*args)
    def erase(*args): return _XmippData.vectori_erase(*args)
    def __init__(self, *args): 
        this = _XmippData.new_vectori(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _XmippData.vectori_push_back(*args)
    def front(*args): return _XmippData.vectori_front(*args)
    def back(*args): return _XmippData.vectori_back(*args)
    def assign(*args): return _XmippData.vectori_assign(*args)
    def resize(*args): return _XmippData.vectori_resize(*args)
    def insert(*args): return _XmippData.vectori_insert(*args)
    def reserve(*args): return _XmippData.vectori_reserve(*args)
    def capacity(*args): return _XmippData.vectori_capacity(*args)
    __swig_destroy__ = _XmippData.delete_vectori
    __del__ = lambda self : None;
vectori_swigregister = _XmippData.vectori_swigregister
vectori_swigregister(vectori)

class vectord(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectord, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectord, name)
    __repr__ = _swig_repr
    def iterator(*args): return _XmippData.vectord_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _XmippData.vectord___nonzero__(*args)
    def __len__(*args): return _XmippData.vectord___len__(*args)
    def pop(*args): return _XmippData.vectord_pop(*args)
    def __getslice__(*args): return _XmippData.vectord___getslice__(*args)
    def __setslice__(*args): return _XmippData.vectord___setslice__(*args)
    def __delslice__(*args): return _XmippData.vectord___delslice__(*args)
    def __delitem__(*args): return _XmippData.vectord___delitem__(*args)
    def __getitem__(*args): return _XmippData.vectord___getitem__(*args)
    def __setitem__(*args): return _XmippData.vectord___setitem__(*args)
    def append(*args): return _XmippData.vectord_append(*args)
    def empty(*args): return _XmippData.vectord_empty(*args)
    def size(*args): return _XmippData.vectord_size(*args)
    def clear(*args): return _XmippData.vectord_clear(*args)
    def swap(*args): return _XmippData.vectord_swap(*args)
    def get_allocator(*args): return _XmippData.vectord_get_allocator(*args)
    def begin(*args): return _XmippData.vectord_begin(*args)
    def end(*args): return _XmippData.vectord_end(*args)
    def rbegin(*args): return _XmippData.vectord_rbegin(*args)
    def rend(*args): return _XmippData.vectord_rend(*args)
    def pop_back(*args): return _XmippData.vectord_pop_back(*args)
    def erase(*args): return _XmippData.vectord_erase(*args)
    def __init__(self, *args): 
        this = _XmippData.new_vectord(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _XmippData.vectord_push_back(*args)
    def front(*args): return _XmippData.vectord_front(*args)
    def back(*args): return _XmippData.vectord_back(*args)
    def assign(*args): return _XmippData.vectord_assign(*args)
    def resize(*args): return _XmippData.vectord_resize(*args)
    def insert(*args): return _XmippData.vectord_insert(*args)
    def reserve(*args): return _XmippData.vectord_reserve(*args)
    def capacity(*args): return _XmippData.vectord_capacity(*args)
    __swig_destroy__ = _XmippData.delete_vectord
    __del__ = lambda self : None;
vectord_swigregister = _XmippData.vectord_swigregister
vectord_swigregister(vectord)

class intP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, intP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, intP, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_intP(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_intP
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.intP_assign(*args)
    def value(*args): return _XmippData.intP_value(*args)
    def cast(*args): return _XmippData.intP_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _XmippData.intP_frompointer
    if _newclass:frompointer = staticmethod(_XmippData.intP_frompointer)
intP_swigregister = _XmippData.intP_swigregister
intP_swigregister(intP)
intP_frompointer = _XmippData.intP_frompointer

class charP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, charP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, charP, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_charP(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_charP
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.charP_assign(*args)
    def value(*args): return _XmippData.charP_value(*args)
    def cast(*args): return _XmippData.charP_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _XmippData.charP_frompointer
    if _newclass:frompointer = staticmethod(_XmippData.charP_frompointer)
charP_swigregister = _XmippData.charP_swigregister
charP_swigregister(charP)
charP_frompointer = _XmippData.charP_frompointer

class doubleP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, doubleP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, doubleP, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_doubleP(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_doubleP
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.doubleP_assign(*args)
    def value(*args): return _XmippData.doubleP_value(*args)
    def cast(*args): return _XmippData.doubleP_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _XmippData.doubleP_frompointer
    if _newclass:frompointer = staticmethod(_XmippData.doubleP_frompointer)
doubleP_swigregister = _XmippData.doubleP_swigregister
doubleP_swigregister(doubleP)
doubleP_frompointer = _XmippData.doubleP_frompointer

class floatP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, floatP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, floatP, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_floatP(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_floatP
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.floatP_assign(*args)
    def value(*args): return _XmippData.floatP_value(*args)
    def cast(*args): return _XmippData.floatP_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _XmippData.floatP_frompointer
    if _newclass:frompointer = staticmethod(_XmippData.floatP_frompointer)
floatP_swigregister = _XmippData.floatP_swigregister
floatP_swigregister(floatP)
floatP_frompointer = _XmippData.floatP_frompointer

class stringP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, stringP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, stringP, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_stringP(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_stringP
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.stringP_assign(*args)
    def value(*args): return _XmippData.stringP_value(*args)
    def cast(*args): return _XmippData.stringP_cast(*args)
    __swig_getmethods__["frompointer"] = lambda x: _XmippData.stringP_frompointer
    if _newclass:frompointer = staticmethod(_XmippData.stringP_frompointer)
stringP_swigregister = _XmippData.stringP_swigregister
stringP_swigregister(stringP)
stringP_frompointer = _XmippData.stringP_frompointer

class string(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, string, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, string, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_string(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_string
    __del__ = lambda self : None;
string_swigregister = _XmippData.string_swigregister
string_swigregister(string)

class Tabsinc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Tabsinc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Tabsinc, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_Tabsinc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_Tabsinc
    __del__ = lambda self : None;
    def __call__(*args): return _XmippData.Tabsinc___call__(*args)
    def filltable(*args): return _XmippData.Tabsinc_filltable(*args)
Tabsinc_swigregister = _XmippData.Tabsinc_swigregister
Tabsinc_swigregister(Tabsinc)

class KaiserBessel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, KaiserBessel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, KaiserBessel, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_KaiserBessel(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_KaiserBessel
    __del__ = lambda self : None;
    def I0table_maxerror(*args): return _XmippData.KaiserBessel_I0table_maxerror(*args)
    def dump_table(*args): return _XmippData.KaiserBessel_dump_table(*args)
    def sinhwin(*args): return _XmippData.KaiserBessel_sinhwin(*args)
    def i0win(*args): return _XmippData.KaiserBessel_i0win(*args)
    def i0win_tab(*args): return _XmippData.KaiserBessel_i0win_tab(*args)
    def get_window_size(*args): return _XmippData.KaiserBessel_get_window_size(*args)
KaiserBessel_swigregister = _XmippData.KaiserBessel_swigregister
KaiserBessel_swigregister(KaiserBessel)

cdf_gauss = _XmippData.cdf_gauss
cdf_tstudent = _XmippData.cdf_tstudent
cdf_FSnedecor = _XmippData.cdf_FSnedecor
icdf_FSnedecor = _XmippData.icdf_FSnedecor
log2 = _XmippData.log2
randomize_random_generator = _XmippData.randomize_random_generator
student_outside_probb = _XmippData.student_outside_probb
student_within_t0 = _XmippData.student_within_t0
student_outside_t0 = _XmippData.student_outside_t0
student_up_to_t0 = _XmippData.student_up_to_t0
student_from_t0 = _XmippData.student_from_t0
chi2_up_to_t0 = _XmippData.chi2_up_to_t0
chi2_from_t0 = _XmippData.chi2_from_t0
rnd_log = _XmippData.rnd_log
class FileName(string):
    __swig_setmethods__ = {}
    for _s in [string]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, FileName, name, value)
    __swig_getmethods__ = {}
    for _s in [string]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, FileName, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_FileName(*args)
        try: self.this.append(this)
        except: self.this = this
    def compose(*args): return _XmippData.FileName_compose(*args)
    def get_root(*args): return _XmippData.FileName_get_root(*args)
    def get_baseName(*args): return _XmippData.FileName_get_baseName(*args)
    def get_number(*args): return _XmippData.FileName_get_number(*args)
    def get_extension(*args): return _XmippData.FileName_get_extension(*args)
    def init_random(*args): return _XmippData.FileName_init_random(*args)
    def add_prefix(*args): return _XmippData.FileName_add_prefix(*args)
    def add_extension(*args): return _XmippData.FileName_add_extension(*args)
    def without_extension(*args): return _XmippData.FileName_without_extension(*args)
    def without_root(*args): return _XmippData.FileName_without_root(*args)
    def insert_before_extension(*args): return _XmippData.FileName_insert_before_extension(*args)
    def remove_extension(*args): return _XmippData.FileName_remove_extension(*args)
    def remove_all_extensions(*args): return _XmippData.FileName_remove_all_extensions(*args)
    def substitute_extension(*args): return _XmippData.FileName_substitute_extension(*args)
    def without(*args): return _XmippData.FileName_without(*args)
    def remove_until_prefix(*args): return _XmippData.FileName_remove_until_prefix(*args)
    def remove_directories(*args): return _XmippData.FileName_remove_directories(*args)
    def __str__(*args): return _XmippData.FileName___str__(*args)
    __swig_destroy__ = _XmippData.delete_FileName
    __del__ = lambda self : None;
FileName_swigregister = _XmippData.FileName_swigregister
FileName_swigregister(FileName)
solve_2nd_degree_eq = _XmippData.solve_2nd_degree_eq
gaussian1D = _XmippData.gaussian1D
tstudent1D = _XmippData.tstudent1D
gaussian2D = _XmippData.gaussian2D
init_random_generator = _XmippData.init_random_generator
rnd_unif = _XmippData.rnd_unif
rnd_student_t = _XmippData.rnd_student_t
rnd_gaus = _XmippData.rnd_gaus
gaus_within_x0 = _XmippData.gaus_within_x0
gaus_outside_x0 = _XmippData.gaus_outside_x0
gaus_up_to_x0 = _XmippData.gaus_up_to_x0
gaus_from_x0 = _XmippData.gaus_from_x0

exists = _XmippData.exists
xmippBaseDir = _XmippData.xmippBaseDir
init_progress_bar = _XmippData.init_progress_bar
progress_bar = _XmippData.progress_bar
TimeMessage = _XmippData.TimeMessage
IsBigEndian = _XmippData.IsBigEndian
IsLittleEndian = _XmippData.IsLittleEndian
divide_equally = _XmippData.divide_equally
class DocLine(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DocLine, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DocLine, name)
    __repr__ = _swig_repr
    NOT_CONSIDERED = _XmippData.DocLine_NOT_CONSIDERED
    NOT_ASSIGNED = _XmippData.DocLine_NOT_ASSIGNED
    DATALINE = _XmippData.DocLine_DATALINE
    COMMENT = _XmippData.DocLine_COMMENT
    def __init__(self, *args): 
        this = _XmippData.new_DocLine(*args)
        try: self.this.append(this)
        except: self.this = this
    def assign(*args): return _XmippData.DocLine_assign(*args)
    def get(*args): return _XmippData.DocLine_get(*args)
    def set(*args): return _XmippData.DocLine_set(*args)
    def get_text(*args): return _XmippData.DocLine_get_text(*args)
    def get_key(*args): return _XmippData.DocLine_get_key(*args)
    def get_no_components(*args): return _XmippData.DocLine_get_no_components(*args)
    def clear(*args): return _XmippData.DocLine_clear(*args)
    def Is_comment(*args): return _XmippData.DocLine_Is_comment(*args)
    def Is_data(*args): return _XmippData.DocLine_Is_data(*args)
    def set_type(*args): return _XmippData.DocLine_set_type(*args)
    def read(*args): return _XmippData.DocLine_read(*args)
    def __str__(*args): return _XmippData.DocLine___str__(*args)
    __swig_destroy__ = _XmippData.delete_DocLine
    __del__ = lambda self : None;
DocLine_swigregister = _XmippData.DocLine_swigregister
DocLine_swigregister(DocLine)
wait_until_stable_size = _XmippData.wait_until_stable_size
create_empty_file = _XmippData.create_empty_file

class DocFile(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DocFile, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DocFile, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_DocFile(*args)
        try: self.this.append(this)
        except: self.this = this
    def clear(*args): return _XmippData.DocFile_clear(*args)
    def reserve(*args): return _XmippData.DocFile_reserve(*args)
    def assign(*args): return _XmippData.DocFile_assign(*args)
    def debug(*args): return _XmippData.DocFile_debug(*args)
    def read(*args): return _XmippData.DocFile_read(*args)
    def append(*args): return _XmippData.DocFile_append(*args)
    def write(*args): return _XmippData.DocFile_write(*args)
    def go_beginning(*args): return _XmippData.DocFile_go_beginning(*args)
    def go_first_data_line(*args): return _XmippData.DocFile_go_first_data_line(*args)
    def adjust_to_data_line(*args): return _XmippData.DocFile_adjust_to_data_line(*args)
    def next(*args): return _XmippData.DocFile_next(*args)
    def previous(*args): return _XmippData.DocFile_previous(*args)
    def next_data_line(*args): return _XmippData.DocFile_next_data_line(*args)
    def jump(*args): return _XmippData.DocFile_jump(*args)
    def search(*args): return _XmippData.DocFile_search(*args)
    def search_comment(*args): return _XmippData.DocFile_search_comment(*args)
    def remove_multiple_strings(*args): return _XmippData.DocFile_remove_multiple_strings(*args)
    def get_selfile(*args): return _XmippData.DocFile_get_selfile(*args)
    def locate(*args): return _XmippData.DocFile_locate(*args)
    def eof(*args): return _XmippData.DocFile_eof(*args)
    def name(*args): return _XmippData.DocFile_name(*args)
    def exists(*args): return _XmippData.DocFile_exists(*args)
    def getColNumberFromHeader(*args): return _XmippData.DocFile_getColNumberFromHeader(*args)
    def FirstLine_colNumber(*args): return _XmippData.DocFile_FirstLine_colNumber(*args)
    def dataLineNo(*args): return _XmippData.DocFile_dataLineNo(*args)
    def LineNo(*args): return _XmippData.DocFile_LineNo(*args)
    def get_last_key(*args): return _XmippData.DocFile_get_last_key(*args)
    def get_current_key(*args): return _XmippData.DocFile_get_current_key(*args)
    def FirstKey(*args): return _XmippData.DocFile_FirstKey(*args)
    def set_FirstKey(*args): return _XmippData.DocFile_set_FirstKey(*args)
    def get_current_valNo(*args): return _XmippData.DocFile_get_current_valNo(*args)
    def __call__(*args): return _XmippData.DocFile___call__(*args)
    def get_angles(*args): return _XmippData.DocFile_get_angles(*args)
    def get_angles1(*args): return _XmippData.DocFile_get_angles1(*args)
    def get_angles2(*args): return _XmippData.DocFile_get_angles2(*args)
    def set_angles(*args): return _XmippData.DocFile_set_angles(*args)
    def get_image(*args): return _XmippData.DocFile_get_image(*args)
    def get_imagename(*args): return _XmippData.DocFile_get_imagename(*args)
    def set(*args): return _XmippData.DocFile_set(*args)
    def get_current_line(*args): return _XmippData.DocFile_get_current_line(*args)
    def renum(*args): return _XmippData.DocFile_renum(*args)
    def remove(*args): return _XmippData.DocFile_remove(*args)
    def remove_current(*args): return _XmippData.DocFile_remove_current(*args)
    def insert_data_line(*args): return _XmippData.DocFile_insert_data_line(*args)
    def insert_comment(*args): return _XmippData.DocFile_insert_comment(*args)
    def insert_line(*args): return _XmippData.DocFile_insert_line(*args)
    def append_data_line(*args): return _XmippData.DocFile_append_data_line(*args)
    def append_angles(*args): return _XmippData.DocFile_append_angles(*args)
    def append_comment(*args): return _XmippData.DocFile_append_comment(*args)
    def append_line(*args): return _XmippData.DocFile_append_line(*args)
    def clean_comments(*args): return _XmippData.DocFile_clean_comments(*args)
    def randomize(*args): return _XmippData.DocFile_randomize(*args)
    def perturb_column(*args): return _XmippData.DocFile_perturb_column(*args)
    def merge(*args): return _XmippData.DocFile_merge(*args)
    def random_discard(*args): return _XmippData.DocFile_random_discard(*args)
    def sort_by_filenames(*args): return _XmippData.DocFile_sort_by_filenames(*args)
    def col(*args): return _XmippData.DocFile_col(*args)
    def row(*args): return _XmippData.DocFile_row(*args)
    def setCol(*args): return _XmippData.DocFile_setCol(*args)
    def __str__(*args): return _XmippData.DocFile___str__(*args)
    __swig_destroy__ = _XmippData.delete_DocFile
    __del__ = lambda self : None;
DocFile_swigregister = _XmippData.DocFile_swigregister
DocFile_swigregister(DocFile)
DOCMERGE_KEEP_OLD = _XmippData.DOCMERGE_KEEP_OLD
DOCMERGE_KEEP_NEW = _XmippData.DOCMERGE_KEEP_NEW
DOCMERGE_SUM_COLUMN = _XmippData.DOCMERGE_SUM_COLUMN
DOCMERGE_ERROR = _XmippData.DOCMERGE_ERROR

read_Euler_document_file = _XmippData.read_Euler_document_file
select_images = _XmippData.select_images
get_subset_docfile = _XmippData.get_subset_docfile
checkAngle = _XmippData.checkAngle
CPPSQLITE_ERROR = _XmippData.CPPSQLITE_ERROR
class CppSQLite3Table(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CppSQLite3Table, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CppSQLite3Table, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_CppSQLite3Table(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_CppSQLite3Table
    __del__ = lambda self : None;
    def assign(*args): return _XmippData.CppSQLite3Table_assign(*args)
    def numFields(*args): return _XmippData.CppSQLite3Table_numFields(*args)
    def numRows(*args): return _XmippData.CppSQLite3Table_numRows(*args)
    def fieldName(*args): return _XmippData.CppSQLite3Table_fieldName(*args)
    def fieldValue(*args): return _XmippData.CppSQLite3Table_fieldValue(*args)
    def getIntField(*args): return _XmippData.CppSQLite3Table_getIntField(*args)
    def getFloatField(*args): return _XmippData.CppSQLite3Table_getFloatField(*args)
    def getStringField(*args): return _XmippData.CppSQLite3Table_getStringField(*args)
    def fieldIsNull(*args): return _XmippData.CppSQLite3Table_fieldIsNull(*args)
    def setRow(*args): return _XmippData.CppSQLite3Table_setRow(*args)
    def finalize(*args): return _XmippData.CppSQLite3Table_finalize(*args)
CppSQLite3Table_swigregister = _XmippData.CppSQLite3Table_swigregister
CppSQLite3Table_swigregister(CppSQLite3Table)

class CppSQLite3DB(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CppSQLite3DB, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CppSQLite3DB, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _XmippData.new_CppSQLite3DB(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_CppSQLite3DB
    __del__ = lambda self : None;
    def open(*args): return _XmippData.CppSQLite3DB_open(*args)
    def close(*args): return _XmippData.CppSQLite3DB_close(*args)
    def tableExists(*args): return _XmippData.CppSQLite3DB_tableExists(*args)
    def execDML(*args): return _XmippData.CppSQLite3DB_execDML(*args)
    def execQuery(*args): return _XmippData.CppSQLite3DB_execQuery(*args)
    def execScalar(*args): return _XmippData.CppSQLite3DB_execScalar(*args)
    def getTable(*args): return _XmippData.CppSQLite3DB_getTable(*args)
    def compileStatement(*args): return _XmippData.CppSQLite3DB_compileStatement(*args)
    def lastRowId(*args): return _XmippData.CppSQLite3DB_lastRowId(*args)
    def interrupt(*args): return _XmippData.CppSQLite3DB_interrupt(*args)
    def setBusyTimeout(*args): return _XmippData.CppSQLite3DB_setBusyTimeout(*args)
    __swig_getmethods__["SQLiteVersion"] = lambda x: _XmippData.CppSQLite3DB_SQLiteVersion
    if _newclass:SQLiteVersion = staticmethod(_XmippData.CppSQLite3DB_SQLiteVersion)
CppSQLite3DB_swigregister = _XmippData.CppSQLite3DB_swigregister
CppSQLite3DB_swigregister(CppSQLite3DB)
CppSQLite3DB_SQLiteVersion = _XmippData.CppSQLite3DB_SQLiteVersion

class MetaData(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MetaData, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MetaData, name)
    __repr__ = _swig_repr
    def setValue(*args): return _XmippData.MetaData_setValue(*args)
    def getValue(*args): return _XmippData.MetaData_getValue(*args)
    def readOldSelFile(*args): return _XmippData.MetaData_readOldSelFile(*args)
    def readOldDocFile(*args): return _XmippData.MetaData_readOldDocFile(*args)
    __swig_setmethods__["activeLabels"] = _XmippData.MetaData_activeLabels_set
    __swig_getmethods__["activeLabels"] = _XmippData.MetaData_activeLabels_get
    if _newclass:activeLabels = _swig_property(_XmippData.MetaData_activeLabels_get, _XmippData.MetaData_activeLabels_set)
    def addObject(*args): return _XmippData.MetaData_addObject(*args)
    def read(*args): return _XmippData.MetaData_read(*args)
    NO_OBJECTS_STORED = _XmippData.MetaData_NO_OBJECTS_STORED
    NO_MORE_OBJECTS = _XmippData.MetaData_NO_MORE_OBJECTS
    NO_OBJECT_FOUND = _XmippData.MetaData_NO_OBJECT_FOUND
    def firstObject(*args): return _XmippData.MetaData_firstObject(*args)
    def nextObject(*args): return _XmippData.MetaData_nextObject(*args)
    def lastObject(*args): return _XmippData.MetaData_lastObject(*args)
    def goToObject(*args): return _XmippData.MetaData_goToObject(*args)
    __swig_destroy__ = _XmippData.delete_MetaData
    __del__ = lambda self : None;
    def __init__(self, *args): 
        this = _XmippData.new_MetaData(*args)
        try: self.this.append(this)
        except: self.this = this
    def assign(*args): return _XmippData.MetaData_assign(*args)
    def write(*args): return _XmippData.MetaData_write(*args)
    def isEmpty(*args): return _XmippData.MetaData_isEmpty(*args)
    def clear(*args): return _XmippData.MetaData_clear(*args)
    def fastSearch(*args): return _XmippData.MetaData_fastSearch(*args)
    def findObjects(*args): return _XmippData.MetaData_findObjects(*args)
    def findObjectsInRange(*args): return _XmippData.MetaData_findObjectsInRange(*args)
    def combine(*args): return _XmippData.MetaData_combine(*args)
    def removeObjects(*args): return _XmippData.MetaData_removeObjects(*args)
    def removeObject(*args): return _XmippData.MetaData_removeObject(*args)
    def setPath(*args): return _XmippData.MetaData_setPath(*args)
    def setComment(*args): return _XmippData.MetaData_setComment(*args)
    def getPath(*args): return _XmippData.MetaData_getPath(*args)
    def getComment(*args): return _XmippData.MetaData_getComment(*args)
    def size(*args): return _XmippData.MetaData_size(*args)
    def getFilename(*args): return _XmippData.MetaData_getFilename(*args)
    def detectObjects(*args): return _XmippData.MetaData_detectObjects(*args)
    def countObjects(*args): return _XmippData.MetaData_countObjects(*args)
MetaData_swigregister = _XmippData.MetaData_swigregister
MetaData_swigregister(MetaData)

get_statistics = _XmippData.get_statistics
MDL_UNDEFINED = _XmippData.MDL_UNDEFINED
MDL_FIRST_LABEL = _XmippData.MDL_FIRST_LABEL
MDL_ANGLEROT = _XmippData.MDL_ANGLEROT
MDL_COMMENT = _XmippData.MDL_COMMENT
MDL_ANGLETILT = _XmippData.MDL_ANGLETILT
MDL_ANGLEPSI = _XmippData.MDL_ANGLEPSI
MDL_IMAGE = _XmippData.MDL_IMAGE
MDL_MICROGRAPH = _XmippData.MDL_MICROGRAPH
MDL_CTFMODEL = _XmippData.MDL_CTFMODEL
MDL_SHIFTX = _XmippData.MDL_SHIFTX
MDL_SHIFTY = _XmippData.MDL_SHIFTY
MDL_SHIFTZ = _XmippData.MDL_SHIFTZ
MDL_ENABLED = _XmippData.MDL_ENABLED
MDL_ORIGINX = _XmippData.MDL_ORIGINX
MDL_ORIGINY = _XmippData.MDL_ORIGINY
MDL_ORIGINZ = _XmippData.MDL_ORIGINZ
MDL_WEIGHT = _XmippData.MDL_WEIGHT
MDL_FLIP = _XmippData.MDL_FLIP
MDL_REF = _XmippData.MDL_REF
MDL_MAXCC = _XmippData.MDL_MAXCC
MDL_SERIE = _XmippData.MDL_SERIE
MDL_PMAX = _XmippData.MDL_PMAX
MDL_CTFINPUTPARAMS = _XmippData.MDL_CTFINPUTPARAMS
MDL_PERIODOGRAM = _XmippData.MDL_PERIODOGRAM
MDL_NMA = _XmippData.MDL_NMA
MDL_LAST_LABEL = _XmippData.MDL_LAST_LABEL
class MetaDataContainer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, MetaDataContainer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, MetaDataContainer, name)
    __repr__ = _swig_repr
    def assign(*args): return _XmippData.MetaDataContainer_assign(*args)
    def __init__(self, *args): 
        this = _XmippData.new_MetaDataContainer(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _XmippData.delete_MetaDataContainer
    __del__ = lambda self : None;
    def addValue(*args): return _XmippData.MetaDataContainer_addValue(*args)
    def getValue(*args): return _XmippData.MetaDataContainer_getValue(*args)
    def valueExists(*args): return _XmippData.MetaDataContainer_valueExists(*args)
    def pairExists(*args): return _XmippData.MetaDataContainer_pairExists(*args)
    def deleteValue(*args): return _XmippData.MetaDataContainer_deleteValue(*args)
    def writeValueToFile(*args): return _XmippData.MetaDataContainer_writeValueToFile(*args)
    def writeValueToString(*args): return _XmippData.MetaDataContainer_writeValueToString(*args)
    __swig_getmethods__["codifyLabel"] = lambda x: _XmippData.MetaDataContainer_codifyLabel
    if _newclass:codifyLabel = staticmethod(_XmippData.MetaDataContainer_codifyLabel)
    __swig_getmethods__["decodeLabel"] = lambda x: _XmippData.MetaDataContainer_decodeLabel
    if _newclass:decodeLabel = staticmethod(_XmippData.MetaDataContainer_decodeLabel)
    __swig_getmethods__["isValidLabel"] = lambda x: _XmippData.MetaDataContainer_isValidLabel
    if _newclass:isValidLabel = staticmethod(_XmippData.MetaDataContainer_isValidLabel)
MetaDataContainer_swigregister = _XmippData.MetaDataContainer_swigregister
MetaDataContainer_swigregister(MetaDataContainer)
MetaDataContainer_codifyLabel = _XmippData.MetaDataContainer_codifyLabel
MetaDataContainer_decodeLabel = _XmippData.MetaDataContainer_decodeLabel
MetaDataContainer_isValidLabel = _XmippData.MetaDataContainer_isValidLabel

class vectorm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vectorm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vectorm, name)
    __repr__ = _swig_repr
    def iterator(*args): return _XmippData.vectorm_iterator(*args)
    def __iter__(self): return self.iterator()
    def __nonzero__(*args): return _XmippData.vectorm___nonzero__(*args)
    def __len__(*args): return _XmippData.vectorm___len__(*args)
    def pop(*args): return _XmippData.vectorm_pop(*args)
    def __getslice__(*args): return _XmippData.vectorm___getslice__(*args)
    def __setslice__(*args): return _XmippData.vectorm___setslice__(*args)
    def __delslice__(*args): return _XmippData.vectorm___delslice__(*args)
    def __delitem__(*args): return _XmippData.vectorm___delitem__(*args)
    def __getitem__(*args): return _XmippData.vectorm___getitem__(*args)
    def __setitem__(*args): return _XmippData.vectorm___setitem__(*args)
    def append(*args): return _XmippData.vectorm_append(*args)
    def empty(*args): return _XmippData.vectorm_empty(*args)
    def size(*args): return _XmippData.vectorm_size(*args)
    def clear(*args): return _XmippData.vectorm_clear(*args)
    def swap(*args): return _XmippData.vectorm_swap(*args)
    def get_allocator(*args): return _XmippData.vectorm_get_allocator(*args)
    def begin(*args): return _XmippData.vectorm_begin(*args)
    def end(*args): return _XmippData.vectorm_end(*args)
    def rbegin(*args): return _XmippData.vectorm_rbegin(*args)
    def rend(*args): return _XmippData.vectorm_rend(*args)
    def pop_back(*args): return _XmippData.vectorm_pop_back(*args)
    def erase(*args): return _XmippData.vectorm_erase(*args)
    def __init__(self, *args): 
        this = _XmippData.new_vectorm(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(*args): return _XmippData.vectorm_push_back(*args)
    def front(*args): return _XmippData.vectorm_front(*args)
    def back(*args): return _XmippData.vectorm_back(*args)
    def assign(*args): return _XmippData.vectorm_assign(*args)
    def resize(*args): return _XmippData.vectorm_resize(*args)
    def insert(*args): return _XmippData.vectorm_insert(*args)
    def reserve(*args): return _XmippData.vectorm_reserve(*args)
    def capacity(*args): return _XmippData.vectorm_capacity(*args)
    __swig_destroy__ = _XmippData.delete_vectorm
    __del__ = lambda self : None;
vectorm_swigregister = _XmippData.vectorm_swigregister
vectorm_swigregister(vectorm)



