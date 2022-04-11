from lsdo_kit.simulation.simulation_model import SimulationModel

class Simulation(object):
    def __init__(self, name, nt=1):
        self.name = name
        self.actuation_dict = {}
        self.solver_dict = {}
        self.mesh_list = []
        self.nt = nt
        self.solver_info = {}

    def add_solver_info(self, name, value):
        print(f'{name}:', type(value))
        self.solver_info[name] = value
    
    def add_solver(self, solver_obj):
        solver_obj_csdl_keys = list(solver_obj.sim_solver_inputs_dict.keys())
        for csdl_input_key in solver_obj_csdl_keys:
            if csdl_input_key not in list(self.solver_info.keys()):
                raise Exception(f'{csdl_input_key} information is missing from the solver info dictionary! Use add_solver_info!!')   
        
        self.solver_dict[solver_obj.name] = solver_obj

        if hasattr(solver_obj, 'mesh'):
            self.mesh_list.append(solver_obj.mesh)
        
    def add_actuations(self, design, actuation_list):
        for actuation in actuation_list:
            # print('NAMES NAMES NAMES: ', actuation.name)
            actuation.axis = design.subtract_pointsets(actuation.pointset2, actuation.origin)
            self.actuation_dict[actuation.name] = actuation

    def set_value(self):
        pass

    def connect(self):
        pass