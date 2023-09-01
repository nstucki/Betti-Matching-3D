import gdb

# Start off with defining the printer as a Python object.

class CubePrinter:

    # The constructor takes the value and stores it for later.

    def __init__(self, val):
        self.val = val

    # The to_string method returns the value of the
    # si_signo attribute of the directory.

    def to_string(self):
        cube = self.val
        x = (int(cube['index']) >> 44) & 0xfffff
        y = (int(cube['index']) >> 24) & 0xfffff
        z = (int(cube['index']) >> 4) & 0xfffff
        type = int(cube['index']) & 0xf
        birth = cube['birth']
        dim = int(cube["dimension"])
        # print(cube['parentCgc'].dereference())
        
        # print(cube['parentCgc'], cube['parentCgc'] != 0)
        if cube['parentCgc'] != 0:
            # get_birth_address = hex(cube['parentCgc'].dereference()['getBirthByCoordinate'].address)
            # print(dir(cube['parentCgc'].dereference()['getBirthByCoordinate'].type))
            # print((cube['parentCgc'].dereference()['getBirthByCoordinate'].type))
            # print((cube['parentCgc'].dereference()['getBirthByCoordinate'].type.pointer()))
            # get_birth_type = cube['parentCgc'].dereference()['getBirthByCoordinate'].type.pointer()
            # get_birth = lambda x, y, z: gdb.parse_and_eval(f"(({get_birth_type})({get_birth_address}))({cube['parentCgc']}, {x}, {y}, {z})")

            parent_cgc_address = hex(cube['parentCgc'])
            parent_cgc_type = cube['parentCgc'].type
            get_birth = lambda x, y, z: gdb.parse_and_eval(f"(({parent_cgc_type})({parent_cgc_address}))->getBirthByCoordinate({x}, {y}, {z})")
        else:
            get_birth = lambda *args, **kwargs: "?"
        # print(dir(cube['parentCgc']))
        # cgc = 
        # dim = (x % 2) + (y % 2) + (z % 2)
        def range_str(index, varies):
          return str(index) if not varies else f"{index}-{index+1}"
        if dim == 0:
          return f"•({get_birth(x,y,z)}) at ({x}, {y}, {z}), birth={birth}"
        if dim == 1:
          return f"━({get_birth(x,y,z)}-{get_birth(x+int(type==0),y+int(type==1),z+int(type==2))}) at ({range_str(x, type==0)}, {range_str(y, type==1)}, {range_str(z, type==2)}), birth={birth}"
        if dim == 2:
          # If type == 0: y, z vary, (x,y,z)-(x,y,z+1);(x,y+1,z)-(x,y+1,z+1)
          # If type == 1: x, z vary, (x,y,z)-(x,y,z+1);(x+1,y,z)-(x+1,y,z+1)
          # If type == 2: x, y vary, (x,y,z)-(x,y+1,z);(x+1,y,z)-(x+1,y+1,z)
          return (f"■({get_birth(x,y,z)}-{get_birth(x,y+int(type==2),z+int(type<2))};{get_birth(x+int(type>0),y+int(type==0),z)}-{get_birth(x+int(type>0),y+int(type!=1),z+int(type<2))})" +
                    f" at ({range_str(x, type!=0)}, {range_str(y, type!=1)}, {range_str(z, type!=2)}), birth={birth}")
        if dim == 3:
          return (f"■³({get_birth(x,y,z)}-{get_birth(x,y,z+1)};{get_birth(x,y+1,z)}-{get_birth(x,y+1,z+1)} ;; {get_birth(x+1,y,z)}-{get_birth(x+1,y,z+1)};{get_birth(x+1,y+1,z)}-{get_birth(x+1,y+1,z+1)})" +
                    f" at ({range_str(x, True)}, {range_str(y, True)}, {range_str(z, True)}), birth={birth}")
        # if type == 0:
        #     return f"•/■³(birth={birth}) at ({x},{y},{z})"
        # if type > 0:
        #     return f"━/■(birth={birth}) at ({x},{y},{z})"

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

    def children (self):
        counter = map (self.format_count, itertools.count())
        # Map over the hash table and flatten the result.
        data = list(self.flatten (map (self.format_one, StdHashtableIterator (self.hashtable()))))
        # data = self.flatten (StdHashtableIterator (self.hashtable()))


        # data = [CubePrinter({'index': d, 'birth': '?'}).to_string() if i % 2 == 0 else str(d) for i, d in enumerate(data)]
        # data = [gdb.Value([], "Cube") if i % 2 == 0 else str(d) for i, d in enumerate(data)]
        # data = [CubePrinter({'index': d, 'birth': '?'}).to_string() if i % 2 == 0 else str(d) for i, d in enumerate(data)]


        # Zip the two iterators together.
        # return data
        display_pairs = [int(data[i+1]) for i in range(0, len(data), 2)]
        counter = [f"{CubePrinter({'index': data[i], 'birth': '?'}).to_string()} [{data[i]}]" for i in range(0, len(data), 2)]
        return zip (counter, display_pairs)
        # return zip(counter, [f"{data[i]}: {data[i+1]}" for i in range(0, len(data), 2)])
        # return zip(counter, ["hallo", "welt"])

    def display_hint (self):
        # return 'map'
        return "array"

# Next, define the lookup function that returns the
# printer object when it receives a siginfo_t.

# The function takes the GDB value-type, which, in
# our example is used to look for the siginfo_t.

class Bla:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        return self.val

def pretty_printer_selector(val):
    if "Cube" in str(val.type): return CubePrinter(val)
    # if str(val.type) == 'int64_t': return CubePrinter({'index': val, 'birth': '?'})
    if "std::unordered_map" in str(val.type) and "unsigned long" in str(val.type):
        return CubeUnorderedMapPrinter(str(val.type), val)

# Finally, append the pretty-printer as object/ function to 
# the list of registered GDB printers.

# pp = gdb.printing.RegexpCollectionPrettyPrinter('BettiMatching3D')
# print(pp.enabled)
# print("##############################")
# gdb.printing.PrettyPrinter()
# pp.add_printer("CubePrinter", "Cube", CubePrinter)
# pp.add_printer("CubeUnorderedMapPrinter", "std::unordered_map", CubeUnorderedMapPrinter)
# pp.add_printer('MyUnorderedMapOfMyClassPrinter', 'std::unordered_map<unsigned long, unsigned long>', CubeUnorderedMapPrinter)
gdb.pretty_printers.insert(0, pretty_printer_selector)

# Our pretty-printer is now available when we debug 
# the inferior program in GDB.