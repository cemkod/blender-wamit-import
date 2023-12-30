bl_info = {
    "name": "WAMIT Importer",
    "blender": (2, 80, 0),
    "version": (1, 0),
    "location": "File > Import > WAMIT GDF (.gdf)",
    "description": "Import WAMIT GDF (.gdf) files",
    "author": "Cem Uzunoğlu",
    "category": "Import-Export",
}

import bpy
from bpy_extras.io_utils import ImportHelper
from math import floor


# 6.1 THE GEOMETRIC DATA FILE
# In the low-order method the wetted surface of a body is represented by an ensemble of
# connected four-sided facets, or panels. The Geometric Data File contains a description
# of this discretized surface, including the body length scale, gravity, symmetry indices, the
# total number of panels specified, and for each panel the Cartesian coordinates x, y, z of
# its four vertices. A panel degenerates to a triangle when the coordinates of two vertices
# coincide. The order in which the panels are defined in the file is unimportant, but each
# panel must be described completely by a set of 12 real numbers (three Cartesian coordinates
# for each vertex) which are listed consecutively, with a line break between the last vertex of
# each panel and the first vertex of the next. The value of gravity serves to define the units
# of length, which apply to the body length scale, panel o↵sets, and to all related parameters
# in the other input files. The coordinate system x, y, z in which the panels are defined is
# referred to as the body coordinate system. The only restrictions on the body coordinate
# system are that it is a right-handed Cartesian system and that the zaxis is vertical and
# positive upward.
# The name of the GDF file can be any legal filename accepted by the operating system,
# with a maximum length of 16 ASCII characters, followed by the extension ‘.gdf’.
# The data in the GDF file can be input in the following form:
# header
# ULEN GRAV
# ISX ISY
# NPAN
# X1(1) Y1(1) Z1(1) X2(1) Y2(1) Z2(1) X3(1) Y3(1) Z3(1) X4(1) Y4(1) Z4(1)
# X1(2) Y1(2) Z1(2) X2(2) Y2(2) Z2(2) X3(2) Y3(2) Z3(2) X4(2) Y4(2) Z4(2)
# .
# .
# .
# . . . . . . . . . . . . X4(NPAN) Y4(NPAN) Z4(NPAN)
# Each line of data indicated above is input by a separate FORTRAN READ statement,
# hence line breaks between data must exist as shown. Additional line breaks between data
# perspective view is from above the free surface, showing portions of the exterior and interior of the cylinder
# (lower and upper portions of the figure, respectively). The view of panel i is from the ‘wet side’, inside
# the fluid domain, so the vertex ordering appears anti-clockwise. The view of panel j is from the ‘dry side’
# outside the fluid domain, so the vertex ordering appears clockwise.
# shown above have no e↵ect on the READ statement, so that for example the user may elect
# to place the twelve successive coordinates for each panel on four separate lines. (However
# the format used above is more ecient regarding storage and access time.)
# Input data must be in the order shown above, with at least one blank space separating
# data on the same line.
# The definitions of each entry in this file are as follows:
# ‘header’ denotes a one-line ASCII header dimensioned CHARACTER⇤72. This line is
# available for the user to insert a brief description of the file, with maximum length 72
# characters.
# ULEN is the dimensional length characterizing the body dimension. This parameter
# corresponds to the quantity L used in Chapter 4 to nondimensionalize the quantities output
# from WAMIT. ULEN can be input in any units of length, meters or feet for example, as
# long as the length scale of all other inputs is in the same units. ULEN must be a positive
# number, greater than 105. An error return and warning statement are generated if the
# last restriction is not satisfied.
# GRAV is the acceleration of gravity, using the same units of length as in ULEN. The units
# of time are always seconds. If lengths are input in meters or feet, input 9.80665 or 32.174,
# respectively, for GRAV.
# ISX, ISY are the geometry symmetry indices which have integer values 0, +1. If ISX
# and/or ISY =1, x = 0 and/or y = 0 is a geometric plane of symmetry, and the input
# data (panel vertex coordinates X,Y,Z and their total number NPAN) are restricted to one
# quadrant or one half of the body, namely the portion x > 0 and/or y > 0. Conversely, if
# ISX=0 and ISY=0, the complete submerged surface of the body must be represented by
# panels.
# ISX = 1: The x = 0 plane is a geometric plane of symmetry.
# ISX = 0: The x = 0 plane is not a geometric plane of symmetry.
# ISY = 1: The y = 0 plane is a geometric plane of symmetry.
# ISY = 0: The y = 0 plane is not a geometric plane of symmetry.
# For all values of ISX and ISY, the (x, y) axes are understood to belong to the body system.
# The panel data are always referenced with respect to this system, even if walls or other
# bodies are present.
# NPAN is equal to the number of panels with coordinates defined in this file, i.e. the
# number required to discretize a quarter, half or the whole of the body surface if there exist
# two, one or no planes of symmetry respectively.
# X1(1), Y1(1), Z1(1) are the (x, y, z) coordinates of vertex 1 of the first panel, X2(1),
# Y2(1), Z2(1) the (x, y, z) coordinates of the vertex 2 of the first panel, and so on. These
# are expressed in the same units as the length ULEN. The vertices must be numbered in
# the counter-clockwise direction when the panel is viewed from the fluid domain, as shown
# in Figure 6.1. The precise format of each coordinate is unimportant, as long as there
# is at least one blank space between coordinates, and the coordinates of the four vertices
# representing a panel are listed sequentially.
# There are two situations when panels lie on the free surface, and thus all four vertices
# are on the free surface: (1) the discretization of a structure which has zero draft over
# part or all of its submerged surface, and (2) the discretization of the interior free surface
# for the irregular frequency removal as described in Chapter 10. For the first case, where
# the panels are part of the physical surface, the panel vertices must be numbered in the
# counter-clockwise direction when the panel is viewed from the fluid domain as in the case
# of submerged panels. For the second case, where the panel is interior to the body and
# non physical, the vertices must be numbered in the clockwise direction when the panel is
# viewed from inside the structure (or in the counter-clockwise direction when the panel is
# viewed from above the free surface). Details of the discretization of the interior free surface
# are provided in Chapter 10.
# Although the panels on the free surface are legitimate in these two special cases, a warning
# message is displayed by WAMIT when it detects panels with zero draft, which have four
# vertices on the free surface. This is to provide a warning to users for a possible error in
# the discretization other than the above two exceptional cases. The run continues in this
# case, without interruption. An error message is displayed with an interruption of the run
# when the panels have only three vertices on the free surface, unless two adjacent vertices
# are coincident. (The latter provision permits the analysis of a triangular panel with one
# side in the free surface.)
# The three Cartesian coordinates of four vertices must always be input for each panel, in
# a sequence of twelve real numbers. Triangles are represented by allowing the coordinates
# of two adjacent vertices to coincide, as in the center bottom panels shown in Figure 6.1.
# Two adjacent vertices are defined to be coincident if their included side has a length less
# than ULEN ⇥ 106. An error return results if the computed area of any panel is less than
# ULEN2 ⇥ 1010.
# The input vertices of a panel do not need to be co-planar. WAMIT internally defines planar
# panels that are a best fit to four vertices not lying on a plane. However it is advisable to
# discretize the body so that the input vertices defining each panel lie close to a plane, in
# order to achieve good accuracy in the computed velocity potentials. An error message is
# printed if a panel has two intersecting sides. A warning message is printed if a panel is
# ‘convex’ (the included angle between two adjacent sides exceeds 180 degrees).
# The origin of the body coordinate system may be on, above or below the free surface. The
# vertical distance of the origin from the free surface is specified in the Potential Control File.
# The same body-system is also used to define the forces, moments, and body motions. (See
# Chapter 5 regarding the change in reference of phase relations when walls are present.)
# Only the wetted surface of the body should be paneled, and then only half or a quarter
# of it if there exist one or two planes of symmetry respectively. This also applies to bodies
# mounted on the sea bottom or on one or two vertical walls. The number of panels NPAN
# refers to the number used to discretize a quarter, half or the whole body wetted surface if
# two, one or no planes of symmetry are present respectively.
# The displaced volume of the structure deserves particular discussion. Three separate algorithms are used to evaluate this quantity, as explained in Section 3.1. Except for the special
# case where the structure is bottom-mounted, the three evaluations (VOLX, VOLY, VOLZ)
# should be identical, but they will generally di↵er by small amounts due to inaccuracies in
# machine computation and, more significantly, to approximations in the discretization of
# the body surface.

# Example file:
# 1 body structure
# 1 9.80665 ULEN GRAV
# 0 0 ISX ISY
#  
# -7.500000 -2.500000 1.000000 -6.500000 -2.500000 1.000000 -6.500000 -1.500000 1.000000 -7.500000 -1.500000 1.000000 
# -6.500000 -2.500000 1.000000 -5.500000 -2.500000 1.000000 -5.500000 -1.500000 1.000000 -6.500000 -1.500000 1.000000 
# -5.500000 -2.500000 1.000000 -4.500000 -2.500000 1.000000 -4.500000 -1.500000 1.000000 -5.500000 -1.500000 1.000000 
# -4.500000 -2.500000 1.000000 -3.500000 -2.500000 1.000000 -3.500000 -1.500000 1.000000 -4.500000 -1.500000 1.000000 
# -3.500000 -2.500000 1.000000 -2.500000 -2.500000 1.000000 -2.500000 -1.500000 1.000000 -3.500000 -1.500000 1.000000 
# -2.500000 -2.500000 1.000000 -1.500000 -2.500000 1.000000 -1.500000 -1.500000 1.000000 -2.500000 -1.500000 1.000000 
# -1.500000 -2.500000 1.000000 -0.500000 -2.500000 1.000000 -0.500000 -1.500000 1.000000 -1.500000 -1.500000 1.000000 


def import_wamit_gdf_operator(context, filepath):
    # print a hello to the blender python console
    print("Hello World")
    # Open the GDF file
    with open(filepath, 'r') as file:
        lines = file.readlines()

    # Initialize default values
    ulen, grav, isx, isy, npan = 0.0, 0.0, 0.0, 0.0, 0.0

    # Create mesh object
    mesh = bpy.data.meshes.new(name='WAMIT_Object')
    obj = bpy.data.objects.new('WAMIT_Object', mesh)

    # Link the object to the scene
    bpy.context.scene.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Create mesh data
    mesh.from_pydata([], [], [])    

    vertices = []
    faces = []

    # Parse header
    # it is header until the first blank line

    keywords = ['ULEN', 'GRAV', 'ISX', 'ISY', 'NPAN']
    in_header = True

    for i in range(len(lines)):
        words = lines[i].strip().split()

        if in_header:

            if len(words) == 1:
                print("End of header")
                in_header = False
                continue

            # loop through keywords
            for keyword in keywords:
                # if keyword is found
                if keyword in words:
                    # get the index of the keyword
                    index = words.index(keyword)
                    valueindex = index - (floor(len(words)/2));
                    # get the value of the keyword
                    value = words[valueindex]
                    # assign the value to the variable
                    if keyword == 'ULEN':
                        print("ULEN: " + value)
                        ulen = float(value)
                    elif keyword == 'GRAV':
                        print("GRAV: " + value)
                        grav = float(value)
                    elif keyword == 'ISX':
                        print("ISX: " + value)
                        isx = float(value)
                    elif keyword == 'ISY':
                        print("ISY: " + value)
                        isy = float(value)
                    elif keyword == 'NPAN':
                        print("NPAN: " + value)
                        npan = float(value)
        else:
            if len(words) == 12:
                print("Vertex Data: " + lines[i])
                vertices.extend([(float(words[j]), float(words[j + 1]), float(words[j + 2])) for j in range(0, 12, 3)])
                faces.append([k for k in range(len(vertices) - 4, len(vertices))])
        
    if isx == 1:
        print("X = 0 is a geometric plane of symmetry")
        #create a mirror modifier
        bpy.ops.object.modifier_add(type='MIRROR')
        bpy.context.object.modifiers["Mirror"].use_axis[0] = True
    
    if isy == 1:
        print("Y = 0 is a geometric plane of symmetry")
        #create a mirror modifier
        bpy.ops.object.modifier_add(type='MIRROR')
        bpy.context.object.modifiers["Mirror"].use_axis[1] = True


    # Add vertices and faces to the mesh
    mesh.from_pydata(vertices, [], faces)
    mesh.update()

    # Set object properties
    obj.scale = (ulen, ulen, ulen)  # Scale the object based on ULEN

    print("WAMIT GDF file imported successfully.")

class IMPORT_WAMIT_GDF_OT_operator(bpy.types.Operator, ImportHelper):
    bl_idname = "import_mesh.wamit_gdf"
    bl_label = "Import WAMIT GDF"
    bl_options = {'REGISTER', 'UNDO'}

    filename_ext = ".gdf"
    filter_glob: bpy.props.StringProperty(default="*.gdf", options={'HIDDEN'})

    def execute(self, context):
        import_wamit_gdf_operator(context, self.filepath)
        return {'FINISHED'}


def menu_func_import(self, context):
    self.layout.operator(IMPORT_WAMIT_GDF_OT_operator.bl_idname, text="WAMIT GDF (.gdf)")


def register():
    bpy.utils.register_class(IMPORT_WAMIT_GDF_OT_operator)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(IMPORT_WAMIT_GDF_OT_operator)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_import)


if __name__ == "__main__":
    register()
