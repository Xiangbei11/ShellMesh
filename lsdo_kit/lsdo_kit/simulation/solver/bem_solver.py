from lsdo_kit.simulation.solver.solver import Solver

class BemSolver(Solver):
    def __init__(self, name, airfoil : str, num_blades : int, num_radial : int, num_tangential : int):
        super().__init__(name)
        self.geometric_shapes_dict = {}
        self.airfoil = airfoil
        self.num_blades = num_blades
        self.num_radial = num_radial
        self.num_tangential = num_tangential
        self.force_vector_origin = []
        self.force_vector_direction = []

        self.sim_solver_inputs_dict = {
            'velocity' : 'velocity',
            'RPM'      : 'RPM_input',
            'altitude' : 'altitude_input',
        }

        self.sim_solver_outputs_dict = {
            'thrust' : 'total_thrust',
            'torque' : 'total_torque'
        }

        self.design_solver_inputs_dict = {} 
    
    def add_geometric_inputs(self, pointset_name, bem_input_name):
        self.design_solver_inputs_dict[pointset_name] = bem_input_name

    def set_geometric_input_shapes(self, name, shape):
        if name in list(self.design_solver_inputs_dict.keys()):
            self.geometric_shapes_dict[name] = shape
        else:
            raise Exception('Viable geometric input names: rotor_radius, chord, or \
            pitch')

    def add_force_vector(self, vector):
        self.force_vector_origin    = vector[0]
        self.force_vector_direction = vector[1]




if __name__ == "__main__":
    vlm = BemSolver('hello')