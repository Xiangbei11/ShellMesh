from lsdo_kit.simulation.solver.solver import Solver

class VlmSolver(Solver):
    def __init__(self, name, mesh=None):
        super().__init__(name)
        self.mesh = mesh
        
        self.sim_solver_inputs_dict = {
            'velocity' : 'all_timestep_frame_vel',
        }

        self.sim_solver_outputs_dict = {
            'lift' : 'L',
            'drag' : 'D',
            'coefficient of lift' : 'CL',
            'coefficient of drag' : 'CD'
        }

if __name__ == "__main__":
    vlm = VlmSolver('hello')

