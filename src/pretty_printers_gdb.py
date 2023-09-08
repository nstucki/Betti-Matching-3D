import gdb

class CubePrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        cube = self.val
        x = (int(cube['index']) >> 44) & 0xfffff
        y = (int(cube['index']) >> 24) & 0xfffff
        z = (int(cube['index']) >> 4) & 0xfffff
        index = cube['index']
        type = int(cube['index']) & 0xf
        birth = f"{float(cube['birth']):.3f}"
        dim = int(cube["dimension"])

        try:
            parent_cgc = gdb.parse_and_eval("&cgc")
        except gdb.error:
            parent_cgc = 0

        if parent_cgc != 0:
            parent_cgc_address = hex(parent_cgc)
            parent_cgc_type = str(parent_cgc.type)
            def get_birth(x, y, z):
                try:
                    return f'{float(gdb.parse_and_eval(f"(({parent_cgc_type})({parent_cgc_address}))")["grid"][x+1][y+1][z+1]):.3f}'
                except gdb.error as e:
                    return "❌"
        else:
            get_birth = lambda *args, **kwargs: "?"

        def range_str(index, varies):
          return str(index) if not varies else f"{index}-{index+1}"
        if dim == 0:
          return (f"•({get_birth(x,y,z)})" +
                    f" at ({x}, {y}, {z}){f', birth={birth}' if str(birth) != '?' else ''}" +
                    f" [{index}]")
        if dim == 1:
          return (f"━({get_birth(x,y,z)} ━ {get_birth(x+int(type==0),y+int(type==1),z+int(type==2))}) at ({range_str(x, type==0)}, {range_str(y, type==1)}, {range_str(z, type==2)})" +
                    f"{f', birth={birth}' if str(birth) != '?' else ''}" +
                    f" [{index}]")
        if dim == 2:
          # If type == 0: y, z vary, (x,y,z)-(x,y,z+1);(x,y+1,z)-(x,y+1,z+1)
          # If type == 1: x, z vary, (x,y,z)-(x,y,z+1);(x+1,y,z)-(x+1,y,z+1)
          # If type == 2: x, y vary, (x,y,z)-(x,y+1,z);(x+1,y,z)-(x+1,y+1,z)
          return (f"■({get_birth(x,y,z)} ━ {get_birth(x,y+int(type==2),z+int(type<2))} ; {get_birth(x+int(type>0),y+int(type==0),z)} ━ {get_birth(x+int(type>0),y+int(type!=1),z+int(type<2))})" +
                    f" at ({range_str(x, type!=0)}, {range_str(y, type!=1)}, {range_str(z, type!=2)}){f', birth={birth}' if str(birth) != '?' else ''}" +
                    f" [{index}]")
        if dim == 3:
          return (f"■³({get_birth(x,y,z)} ━ {get_birth(x,y,z+1)} ; {get_birth(x,y+1,z)} ━ {get_birth(x,y+1,z+1)} ;; {get_birth(x+1,y,z)} ━ {get_birth(x+1,y,z+1)} ; {get_birth(x+1,y+1,z)} ━ {get_birth(x+1,y+1,z+1)})" +
                    f" at ({range_str(x, True)}, {range_str(y, True)}, {range_str(z, True)}){f', birth={birth}' if str(birth) != '?' else ''}" +
                    f" [{index}]")
        if dim == -1:
            return f"•/━/■/■³(?) at ({x}, {y}, {z}), birth={birth}"

_versioned_namespace = '__8::'
def strip_versioned_namespace(typename):
    global _versioned_namespace
    if _versioned_namespace:
        return typename.replace(_versioned_namespace, '')
    return typename

import itertools

def lookup_templ_spec(templ, *args):
    """
    Lookup template specialization templ<args...>
    """
    t = '{}<{}>'.format(templ, ', '.join([str(a) for a in args]))
    try:
        return gdb.lookup_type(t)
    except gdb.error as e:
        # Type not found, try again in versioned namespace.
        global _versioned_namespace
        if _versioned_namespace and _versioned_namespace not in templ:
            t = t.replace('::', '::' + _versioned_namespace, 1)
            try:
                return gdb.lookup_type(t)
            except gdb.error:
                # If that also fails, rethrow the original exception
                pass
        raise e

class StdHashtableIterator(object):
    def __init__(self, hashtable):
        self.node = hashtable['_M_before_begin']['_M_nxt']
        valtype = hashtable.type.template_argument(1)
        cached = hashtable.type.template_argument(9).template_argument(0)
        node_type = lookup_templ_spec('std::__detail::_Hash_node', str(valtype),
                                      'true' if cached else 'false')
        self.node_type = node_type.pointer()

    def __iter__(self):
        return self

    def __next__(self):
        if self.node == 0:
            raise StopIteration
        elt = self.node.cast(self.node_type).dereference()
        self.node = elt['_M_nxt']
        valptr = elt['_M_storage'].address
        valptr = valptr.cast(elt.type.template_argument(0).pointer())
        return valptr.dereference()


disabled_std_unordered_map_printer = False

class CubeUnorderedMapPrinter:
    "Print a std::unordered_map or tr1::unordered_map"

    def __init__ (self, typename, val):
        self.typename = strip_versioned_namespace(typename)
        self.val = val

    def hashtable (self):
        if self.typename.startswith('std::tr1'):
            return self.val
        return self.val['_M_h']

    def to_string (self):
        count = self.hashtable()['_M_element_count']
        return self.typename

    @staticmethod
    def flatten (list):
        for elt in list:
            for i in elt:
                yield i

    @staticmethod
    def format_one (elt):
        return (elt['first'], elt['second'])

    @staticmethod
    def format_count (i):
        return '[%d]' % i

    # @staticmethod
    # def get_birth(parent_cgc)

    def children (self):
        counter = map (self.format_count, itertools.count())
        # Map over the hash table and flatten the result.
        data = list(self.flatten (map (self.format_one, StdHashtableIterator (self.hashtable()))))
        display_pairs = [int(data[i+1]) for i in range(0, len(data), 2)]

        try:
            parent_cgc = gdb.parse_and_eval("&cgc")
        except gdb.error:
            parent_cgc = 0

        counter = [f"{CubePrinter({'index': data[i], 'birth': '?', 'dimension': 1, 'parentCgc': parent_cgc}).to_string()} [{data[i]}]" for i in range(0, len(data), 2)]
        
        return zip (counter, display_pairs)
        # return zip(counter, [f"{data[i]}: {data[i+1]}" for i in range(0, len(data), 2)])

    def display_hint(self):
        # return 'map'
        return "array"


def pretty_printer_selector(val):
    global disabled_std_unordered_map_printer
    if not disabled_std_unordered_map_printer:
        disabled_std_unordered_map_printer = True
        gdb.execute("disable pretty-printer .*libstdc++.* libstdc++.*;std..unordered_map")
    if "Cube" in str(val.type): return CubePrinter(val)
    # if str(val.type) == 'int64_t': return CubePrinter({'index': val, 'birth': '?'})
    if "std::unordered_map" in str(val.type) and "unsigned long" in str(val.type):
        return CubeUnorderedMapPrinter(str(val.type), val)

gdb.pretty_printers.insert(0, pretty_printer_selector)
