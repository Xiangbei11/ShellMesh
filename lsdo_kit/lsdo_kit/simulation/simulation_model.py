import csdl
import numpy as np
import warnings
# from VLM_package.VLM_system.vlm_system import VLMSystemModel
from lsdo_kit.design.design_geometry.actuation_model import PrescribedActuationModel
from lsdo_kit.simulation.solver.vlm_solver import VlmSolver
from lsdo_kit.simulation.solver.bem_solver import BemSolver
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets
from lsdo_kit.design.design_geometry.mesh_evaluation_model import MeshEvaluationModel
from lsdo_kit.simulation.solver.vlm_solver_model import VLMSolverModel
from lsdo_kit.simulation.solver.bem_solver_model import BEMSolverModel
# from lsdo_kit.design.design_geometry.actuation_model_temp import ActuationModelTemp
# from lsdo_kit.simulation.solver.solver_model import SolverModel



class SimulationModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('simulation_name')
        self.parameters.declare('simulation_obj')
        self.parameters.declare('design')

    def define(self):
        simulation_name = self.parameters['simulation_name']
        simulation_obj = self.parameters['simulation_obj']
        design = self.parameters['design']
        nt = simulation_obj.nt

        ''' 
        Sets up all of the Operating conditions before the solvers and actuation models are called 
        This prevents unwanted feedback from occuring in the model
        '''
        for name, value in simulation_obj.solver_info.items():
            if isinstance(value, np.ndarray):
                if value.shape == (3,):
                    value = np.tile(value, (nt,1))
                    self.create_input(name, val=value, shape=(nt, 3))

                elif value.shape == (nt,3): 
                    self.create_input(name, val=value, shape=(nt, 3))

                else:
                    raise Exception('Velocity must be supplied for every time step')

            elif type(value) is str:
                if name == 'velocity':
                    self.declare_variable(name, shape=(nt,3))
                else:
                    self.declare_variable(name, shape=(nt,))

            else:
                self.create_input(name, val=value, shape=(nt,))
        
        '''
        Set up the actuation model defined for the simulation. 
        Every simulation should have an actuation defined
        '''
        if simulation_obj.actuation_dict:
                actuation_model = PrescribedActuationModel(design=design, actuation_dict=simulation_obj.actuation_dict, nt=nt)
                self.add(actuation_model,'actuation_model', promotes=[])

                self.add(MeshEvaluationModel(design=design, meshes=simulation_obj.mesh_list, nt=nt), 'mesh_eval_model', promotes=[])
                self.connect('actuation_model.actuated_control_points', 'mesh_eval_model.control_points')       
        else:
            warnings.warn('No actuations defined for this simulation!')

        '''
        Check if there are solvers added to the simulation
        Then loop over the solver dictionary 
        Setting it up this way makes sure that the solvers are added last in the model
        '''
        if simulation_obj.solver_dict:
            solver_dictionary = simulation_obj.solver_dict
            for solver_name, solver_obj in solver_dictionary.items():
                pass
                if type(solver_obj) is VlmSolver:
                    # TODO: Discuss if we want vlmsolver to be able to store more than 1 mesh (surface)
                    mesh_name = solver_obj.mesh.name
                    mesh_shape = solver_obj.mesh.pointset_list[0].shape
                    pointset_name = solver_obj.mesh.pointset_list[0].name

                    vlm_model = VLMSolverModel(surface_names=[mesh_name], surface_shapes=[mesh_shape], nt=nt)
                    self.add(vlm_model, solver_name, promotes=[])
                    self.connect(f'mesh_eval_model.{mesh_name}_{pointset_name}', f'{solver_name}.all_timestep_'+mesh_name)  
                
                elif type(solver_obj) is BemSolver:

                    bem_model = BEMSolverModel(airfoil=solver_obj.airfoil,
                                            geometric_shapes_dict = solver_obj.geometric_shapes_dict,
                                            nt = nt,
                                            num_radial = solver_obj.num_radial,
                                            num_tangential = solver_obj.num_tangential,
                                            num_blades = solver_obj.num_blades)
                                            
                    self.add(bem_model, solver_name, promotes=[])
    

        for solver_info_name in simulation_obj.solver_info.keys():
            for solver_name, solver_obj in simulation_obj.solver_dict.items():
                for csdl_input_key, csdl_input_connection in solver_obj.sim_solver_inputs_dict.items():
                    if solver_info_name == csdl_input_key:
                        print('SOLVER INFO NAME: ', solver_info_name)
                        self.connect(solver_info_name, solver_name + '.' + csdl_input_connection)
                    else:
                        continue

        # '''
        # Issue the connections between operating conditions and solvers after solvers have been added
        # '''
        # for solver_name, solver_obj in solver_dictionary.items():
        #     if type(solver_obj) is VlmSolver:
        #         for solver_oc_name, solver_oc_connect in solver_obj.operating_conditions_dict.items():
        #             # for t in range(nt):
        #             self.connect(solver_oc_name, f'all_timestep_' + solver_oc_connect)

        
        # self.connect('mesh_eval_model.vlm_mesh_camber_surface_mesh', 'vlm_model.wing')

        





        # for attribute, value in simulation_obj.__dict__.items():
        #     print('ATTRIBUTES: ', attribute)
        #     if attribute == 'actuation_dict':
        #         actuation_model = PrescribedActuationModel(design=design, actuation_dict=simulation_obj.actuation_dict, nt=simulation_obj.nt)
        #         self.add(actuation_model,'actuation_model', promotes=[])

        #         self.add(MeshEvaluationModel(design=design, meshes=simulation_obj.mesh_list, nt=simulation_obj.nt), 'mesh_eval_model', promotes=[])
        #         self.connect('actuation_model.actuated_control_points', 'mesh_eval_model.control_points')

        #     elif attribute == 'solver_dict':
        #         solver_dictionary = simulation_obj.solver_dict
        #         for solver_name, solver_obj in solver_dictionary.items():
        #             if type(solver_obj) is VlmSolver:

        #                 print('MESH SHAPE: ', solver_obj.mesh.pointset_list[0].shape)
        #                 vlm_model = VLMSolverModel(surface_names=[solver_obj.mesh.name], surface_shapes=[(3,20,3)])
        #                 self.add(vlm_model, 'vlm_model', promotes=[])
        #                 # self.connect('mesh_eval_model.vlm_mesh_camber_surface_mesh', 'vlm_model.wing')
        #                 self.connect('mesh_eval_model.vlm_mesh_camber_surface_mesh', 'vlm_model.'+ solver_obj.mesh.name)

            # elif type(value) is str:
            #     self.declare_variable(simulation_name + '_' + attribute, shape=(1,))
            #     # The value associte with this attribute is going to be used to make a connection 
            #     # NOTE: IF USER IS PROVIDING A CONNECTION THE SHAPE NEEDS TO BE SPECIFIED. WORK IN PROGRESS.

            # elif isinstance(value, numbers.Number):
            #     print('Velocity: ', value)
            #     self.create_input(simulation_name + '_' + attribute, val=np.array([value,0,1]))

            # else:
            #     pass
                # self.create_input(simulation_name + '_' + attribute, val=value)


        # self.connect(simulation_name + '_' + 'velocity', 'vlm_model.frame_vel')

