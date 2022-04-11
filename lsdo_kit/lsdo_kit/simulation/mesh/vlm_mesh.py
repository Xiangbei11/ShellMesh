from lsdo_kit.simulation.mesh.mesh import Mesh

class VLMMesh(Mesh):
    def __init__(self, name, pointset_list=[]) -> None:
        super().__init__(name, pointset_list)