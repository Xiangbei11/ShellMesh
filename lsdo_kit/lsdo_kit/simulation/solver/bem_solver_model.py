from lsdo_rotor.core.BEM_group import BEMGroup
from lsdo_rotor.airfoil.get_surrogate_model import get_surrogate_model
from lsdo_rotor.functions.get_rotor_dictionary import get_rotor_dictionary


from lsdo_rotor.inputs.external_inputs_group import ExternalInputsGroup
from lsdo_rotor.inputs.core_inputs_group import CoreInputsGroup
from lsdo_rotor.inputs.preprocess_group import PreprocessGroup
from lsdo_rotor.core.atmosphere_group import AtmosphereGroup
from lsdo_rotor.core.phi_bracketed_search_group import PhiBracketedSearchGroup
from lsdo_rotor.core.induced_velocity_group import InducedVelocityGroup

from lsdo_rotor import *

from csdl import Model
import csdl
import numpy as np

class BEMSolverModel(Model):
    def initialize(self):
        self.parameters.declare('airfoil')
        self.parameters.declare('geometric_shapes_dict')
        self.parameters.declare('nt')
        self.parameters.declare('num_radial')
        self.parameters.declare('num_tangential')
        self.parameters.declare('num_blades')

    def define(self):
        airfoil = self.parameters['airfoil']
        geometric_shape_dict = self.parameters['geometric_shapes_dict']
        nt = self.parameters['nt']
        num_radial = self.parameters['num_radial']
        num_tangential = self.parameters['num_tangential']
        num_blades = self.parameters['num_blades']
        shape = (nt, num_radial, num_tangential)

        reference_chord     = 0.15
        reference_twist     = 40
        reference_radius    = 0.652

        root_chord = 0.2 * 1.117092 
        tip_chord  = 0.083572 * 1.117092 

        # root_chord = 0.1 
        # tip_chord  = 0.1

        chord = np.linspace(
            tip_chord,
            root_chord,
            num_radial,
        )

        self.create_input('chord', chord)

        root_twist = 41.627000
        tip_twist  = 14.481000
        
        # root_twist          = 80   
        # tip_twist           = 20  

        pitch = np.linspace(
            root_twist * np.pi / 180.,
            tip_twist * np.pi / 180.,
            num_radial,
        )

        self.create_input('pitch', pitch)


        # airfoil = airfoil # 'Clark_Y'
        interp  = get_surrogate_model(airfoil)

        altitude = self.declare_variable('altitude_input', shape=(nt,))
        altitude = csdl.expand(altitude, shape, 'i->ijk')
        
        self.register_output('altitude', altitude)

        V = self.declare_variable('velocity', shape=(nt,3))

        RPM = self.declare_variable('RPM_input', shape=(nt,))
        print(RPM.shape)

        rotor_radius = self.declare_variable('rotor_radius', val=1.105)

        rotational_speed = RPM/60
        hub_radius = 0.2 * rotor_radius
        slice_thickness = (rotor_radius - hub_radius) / (num_radial - 1)

        self.register_output('hub_radius', hub_radius)
        self.register_output('rotational_speed', rotational_speed)
        self.register_output('slice_thickness', slice_thickness)

        reference_rotational_speed = RPM/60

        self.create_input('reference_radius', reference_radius)
        self.create_input('reference_position', shape=(1,3))
        self.create_input('reference_inflow_velocity', val=np.ones((1,1,1,3)))
        self.create_input('reference_pitch', shape=1)
        self.create_input('reference_chord', reference_chord)
        self.create_input('reference_twist', reference_twist)
        self.create_input('reference_axial_inflow_velocity', shape=1)
        self.create_input('reference_blade_solidity', num_blades * reference_chord / 2 / np.pi / reference_radius)
        self.create_input('reference_tangential_inflow_velocity', shape=1)
        self.register_output('reference_rotational_speed', reference_rotational_speed)
        

        self.create_input('position', shape=(nt, 3))
        self.create_input('x_dir', val=np.tile([1,0,0],(nt,1)), shape=(nt, 3))
        self.create_input('y_dir', val=np.tile([0,1,0],(nt,1)), shape=(nt, 3))
        self.create_input('z_dir', val=np.tile([0,0,1],(nt,1)), shape=(nt, 3))



        Vx = V[:,0]
        Vx = csdl.expand(Vx, (nt, num_radial, num_tangential,1), 'ij->iklj')
        Vy = V[:,1]
        Vy = csdl.expand(Vy, (nt, num_radial, num_tangential,1), 'ij->iklj')
        Vz = V[:,2]
        Vz = csdl.expand(Vz, (nt, num_radial, num_tangential,1), 'ij->iklj')
        reference_inflow_velocity = self.create_output('inflow_velocity', shape=shape + (3,))
        reference_inflow_velocity[:,:,:,0] = Vx
        reference_inflow_velocity[:,:,:,1] = Vy
        reference_inflow_velocity[:,:,:,2] = Vz
    

        self.add(ExternalInputsGroup(
            shape = shape,
            num_evaluations = nt,
            num_radial = num_radial,
            num_tangential = num_tangential,
        ), name = 'external_inputs_group')#, promotes = ['*'])
    
        self.add(CoreInputsGroup(
            num_evaluations=nt,
            num_radial=num_radial,
            num_tangential=num_tangential,
        ), name = 'core_inputs_group')#, promotes=['*'])

        self.add(PreprocessGroup(
            num_blades = num_blades,
            shape = shape,
        ), name = 'preprocess_group')#, promotes = ['*'])

        self.add(AtmosphereGroup(
            shape = shape,
        ), name = 'atmosphere_group', promotes = ['*'])

        self.add(PhiBracketedSearchGroup(
            interp = interp,
            shape = shape,
            num_blades = num_blades,
            mode = 2,
        ), name = 'phi_bracketed_search_group')#, promotes = ['*'])

        self.add(InducedVelocityGroup(
            num_blades = num_blades,
            mode  = 2,
            shape = shape,
        ), name = 'induced_velocity_group')#, promotes = ['*'])
