from locale import normalize
import csdl
from csdl_om import Simulator
import numpy as np

from lsdo_kit.design.design_geometry.core.actuation import Actuation
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets
from lsdo_kit.quat_rotate_vector_custom import QuatRotateVectorCustom
from scipy.sparse import coo_array

class PrescribedActuationModel(csdl.Model):
    '''
    Maps from the design geometry to the actuated geometry meshes.
    # Maybe PrescribedActuationModel or other ideas

    Notes:
    -Expands meshes to have a time dimension (nx,ny,3 --> nt,nx,ny,3)
    -Also generates a velocity mesh (nt,nx,ny,3)? Actuation does cause velocity
    -Will only allow bulk actuation
    -We are ok with a parent/child relationship
    --It creates a graph, and evluates in a topological order
    --We use pointsets to define the origin and axis (we subtract 2 pointsets to get the axis)

    -We need to develop API to specify whether actuation profile will be DV vs. input.
    
    
    # Warning: I forget what I was implying with this.
    -Do we want to make the input design a constant input that must be fed in?
      -I can't imagine you would want to make all of the cp (x,y,z) design vars

    '''

    def initialize(self):
      self.parameters.declare('design')
      # self.parameters.declare('meshes')
      self.parameters.declare('actuation_dict')
      self.parameters.declare('nt')


    def define(self):
      design = self.parameters['design']
      nt = self.parameters['nt']
      actuation_dict = self.parameters['actuation_dict']

      geo = design.design_geometry

      pre_actuated_control_points = self.create_input('total_geometry_control_points', val=geo.total_cntrl_pts_vector)
      pre_actuated_control_points_expanded = csdl.expand(pre_actuated_control_points, (nt,) + pre_actuated_control_points.shape, 'ij->tij')

      ''' Copying and Pasting the Evaluation in here  '''
      eval_map = self.create_input(f'eval_map', val=design.design_geometry.eval_map.todense())
      offset_map = self.create_input(f'offset_map', val=design.design_geometry.offset_eval_map)

      num_points = pre_actuated_control_points.shape[0]

      
      # current_actuated_control_points = pre_actuated_control_points_expanded[t,:,:]

      profile_list = []
      # topological order must be enforced
      for actuation_name, actuation_obj in actuation_dict.items():
          profile = actuation_obj.actuation_profile
          # print('profile: ', type(profile))
          if profile is type(str):
            profile_list.append(self.declare_variable(profile, shape=(nt,)))

          elif type(profile) is np.ndarray:
            if profile.shape == (nt,) or profile.shape == (1,):
              profile_list.append(self.create_input(actuation_name, shape=(nt,), val=profile))
            else:
              raise Exception('Actuation profile must be of length nt')
          else:
            raise Exception('Actuation profile must be a str to specify a connection, or a numpy array.')

      act_counter = 0
      for actuation_name, actuation_obj in actuation_dict.items():
        
        actuating_output_indices = None
        for component in actuation_obj.actuating_components:
          for embedded_entity in component.embedded_entities_objects:

            actuating_entity_indices = np.arange(np.cumprod(embedded_entity.shape)[-2]) + embedded_entity.starting_geometry_index

            if actuating_output_indices is None:
              actuating_output_indices = actuating_entity_indices
            else:
              actuating_output_indices = np.vstack((actuating_output_indices, actuating_entity_indices))

        actuating_output_indices = actuating_output_indices.flatten()
        num_ind = len(actuating_output_indices)

        if num_ind == num_points:
          data = np.ones((num_ind))
          sprs = coo_array((data, (np.arange(num_ind), actuating_output_indices)), shape=(num_ind, num_points))


        else:
          # print('IS THIS EVEN RUNNING?!?!?')
          data = np.ones((num_ind))
          sprs = coo_array((data, (np.arange(num_ind), actuating_output_indices)), shape=(num_ind, num_points))


          non_actuating_output_indices = np.delete(np.arange(num_points), actuating_output_indices)
          num_ind_non_act = len(non_actuating_output_indices)
          data = np.ones((num_ind_non_act))
          sprs_mat_non_act = coo_array((data, (np.arange(num_ind_non_act), non_actuating_output_indices)), shape=(num_ind_non_act, num_points))


        all_points = csdl.matmat(eval_map, pre_actuated_control_points)

        origin = actuation_obj.origin
        axis = actuation_obj.axis
        axis_point = all_points[int(axis.output_starting_ind), :]
        origin_point = all_points[int(origin.output_starting_ind), :]

        normalize_axis = axis_point / csdl.expand(csdl.pnorm(axis_point), axis_point.shape, 'i->ij')
        
        thetas = profile_list[act_counter]

        quat = self.create_output(f'quat_{act_counter}', shape= (nt, num_ind) + (4,))
          

        # print('NORMALIZE AXIS SHAPE: ', normalize_axis[0,0].shape)
        # print('THETAS shape: ', thetas.shape)    
    
        normalize_axis_x = csdl.expand( csdl.reshape(normalize_axis[0,0], new_shape=(1,)) , thetas.shape, 'i->j')
        normalize_axis_y = csdl.expand( csdl.reshape(normalize_axis[0,1], new_shape=(1,)) , thetas.shape, 'i->j')
        normalize_axis_z = csdl.expand( csdl.reshape(normalize_axis[0,2], new_shape=(1,)) , thetas.shape, 'i->j')

        
        quat[:,:,0] = csdl.expand(csdl.cos(thetas/2), (nt, num_ind) + (1,), 't->tij')
        quat[:,:,1] = csdl.expand(csdl.sin(thetas/2) * normalize_axis_x, (nt, num_ind) + (1,), 't->tij')
        quat[:,:,2] = csdl.expand(csdl.sin(thetas/2) * normalize_axis_y, (nt, num_ind) + (1,), 't->tij')
        quat[:,:,3] = csdl.expand(csdl.sin(thetas/2) * normalize_axis_z, (nt, num_ind) + (1,), 't->tij')

        actuated_control_points = self.create_output(f'actuated_control_points', val= np.tile(geo.total_cntrl_pts_vector, (nt,1,1)))

        for t in range(nt):
   
          if act_counter == 0 and num_points == num_ind:
            # print('ACT COUNTER 0, EQUAL INDS AND POINTS')
            actuating_points = csdl.sparsematmat(csdl.reshape(pre_actuated_control_points_expanded[t,:,:], new_shape=(num_points, 3)), sparse_mat=sprs )
          
          elif act_counter == 0 and num_points != num_ind:
            # print('ACT COUNTER 0, DIFFERENT INDS AND POINTS')
            # print('SPARSE MATRIX SHAPE: ', sprs.shape)
            # print('PREACTUATED CONTROL POINTS SHAPE: ', pre_actuated_control_points_expanded[t,:,:].shape)
            temp3 = csdl.reshape(pre_actuated_control_points_expanded[t,:,:], new_shape=(num_points, 3))

            self.register_output(f'temp3_{t}', temp3)
            
            actuating_points = csdl.sparsematmat(temp3, sparse_mat=sprs )
            # print('ACTUATING POINTS SHAPE IN IF: ', actuating_points.shape)

            temp3 = csdl.reshape(pre_actuated_control_points_expanded[t,:,:], new_shape=(num_points, 3))
            non_actuating_points = csdl.sparsematmat(temp3, sparse_mat=sprs_mat_non_act )
            self.register_output(f'non_actuating_pt_{t}', non_actuating_points)
            # print('NON SPARSE MATRIX SHAPE: ', sprs_mat_non_act.shape)
            # print('NON ACTUATING POINTS SHAPE: ', non_actuating_points.shape)

          else:
            actuating_points = csdl.sparsematmat(actuated_control_points[t,:,:], sparse_mat=sprs)
            non_actuating_points = csdl.sparsematmat(actuated_control_points[t,:,:], sparse_mat=sprs_mat_non_act)

          if t == 0:
            origin_point = csdl.reshape(origin_point, new_shape=(3,))
            origin_point = csdl.expand(origin_point, actuating_points.shape, 'i->ji')

          self.register_output(f'actuating_pt_{t}', actuating_points)
          
          control_points_origin_frame = actuating_points - origin_point

          self.register_output(f'cp_origin_frame_{t}', control_points_origin_frame)


          if actuation_obj.rot_or_trans == 'rot':
            # print()
            # print('TIMESTEP: ', t)
            # print('NUMBER OF INDICES: ', num_ind)
            # print('NUMBER OF QUATERNIONS: ', quat.shape)
            # print('ACTUATING POINTS SHAPE: ', actuating_points.shape)
            # print('CONTROL POINTS SHAPE: ', control_points_origin_frame.shape)

            temp_quat = csdl.reshape(quat[t,:,:], new_shape=(num_ind, 4))
            rotated_control_points_origin_frame = csdl.quatrotvec(temp_quat, control_points_origin_frame)
            origin_point = csdl.reshape(origin_point, new_shape=rotated_control_points_origin_frame.shape)
            
            temp_control_points = rotated_control_points_origin_frame + origin_point
            self.register_output(f'temp_cp_{t}', temp_control_points)

            if act_counter == 0 and num_points == num_ind:
              actuated_control_points[t,:,:] = csdl.reshape( csdl.sparsematmat(temp_control_points, sparse_mat=sprs.transpose()) , new_shape=(1,num_points,3))

            elif act_counter == 0 and num_points != num_ind:

              actuated_control_points[t,:,:] = csdl.reshape( csdl.sparsematmat(temp_control_points, sparse_mat=sprs.transpose()) + \
                csdl.sparsematmat(non_actuating_points, sparse_mat=sprs_mat_non_act.transpose()) , new_shape = (1, num_points, 3))



          # actuating_points = self.create_output(f'{actuation_name}_points', shape=(nt,) + (num_ind,3))
          # for i in range(3):
          #   actuating_points[:,:,i] = csdl.transpose(csdl.matmat(sprs_mat, csdl.transpose(pre_actuated_control_points_expanded[:,:,i])))



            # pre_actuated_control_points_expanded = 
            # self.register_output('actuated_control_points', actuated_control_points)
            # act_counter += 1























      # profile_list = []
      # val_list = []
      # print('STARTING THE ACTUATION MODEL')
      # # topological order must be enforced
      # for actuation_name, actuation_obj in actuation_dict.items():
      #     profile = actuation_obj.actuation_profile
      #     # print('profile: ', type(profile))
      #     if profile is type(str):
      #       profile_list.append(self.declare_variable(profile, shape=(nt,)))

      #     elif type(profile) is np.ndarray:
      #       if profile.shape == (nt,) or profile.shape == (1,):
      #         profile_list.append(self.create_input(actuation_name, shape=(nt,), val=profile))
      #       else:
      #         raise Exception('Actuation profile must be of length nt')
      #     else:
      #       raise Exception('Actuation profile must be a str to specify a connection, or a numpy array.')




      # '''
      # # TODO: CONFIRM THAT CONTROL POINTS COMING FROM EMBEDDED ENTITIES APPEAR IN SEQUENCIAL ORDER IN ORDER TO USE SLICES
      # '''
      # num_of_actuations = len(actuation_dict.keys())
      # act_ind_dict = {}
      # comp_ind_list = []

      # for actuation_name, actuation_obj in actuation_dict.items():
      #   for component in actuation_obj.actuating_components:
      #     actuating_output_indices = None
      #     for embedded_entity in component.embedded_entities_objects:
      #       actuating_entity_indices = np.arange(np.cumprod(embedded_entity.shape)[1]) + embedded_entity.starting_geometry_index

      #       if actuating_output_indices is None:
      #         actuating_output_indices = actuating_entity_indices
      #       else:
      #         actuating_output_indices = np.vstack((actuating_output_indices, actuating_entity_indices))
          
      #     comp_ind_list.append((np.min(actuating_output_indices), np.max(actuating_output_indices)))
      #   act_ind_dict[actuation_name] = (comp_ind_list)

      # # for t in np.arange(2).astype(int):
      # for t in range(nt):  # This is what originally was here, I replaced it cause I think nt is set up incorrectly rn
      #   act_counter = 0
      #   current_actuated_control_points = pre_actuated_control_points_expanded[t,:,:]
      #   outer_act_list = []
      #   for actuation_name, actuation_obj in actuation_dict.items():

      #     if actuation_obj is not type(Actuation):
      #       Exception('Please pass in a list of actuations.')

      #     origin = actuation_obj.origin
      #     axis = actuation_obj.axis
            
      #     current_actuated_control_points = csdl.reshape(current_actuated_control_points, new_shape=(num_points, 3))
      #     all_points = csdl.matmat(eval_map, current_actuated_control_points)
        

      #     thetas = profile_list[act_counter]

      #     ''' 
      #     We pull out the indices corresponding to the component that needs to be actuated. 
      #     We do this so we can slice the indices
      #     '''
      #     # for component in actuation_obj.actuating_components:
      #     #   for embedded_entity in component.embedded_entities_objects:

      #     #     actuating_entity_indices = np.arange(np.cumprod(embedded_entity.shape)[1]) + embedded_entity.starting_geometry_index

      #     #     if actuating_output_indices is None:
      #     #       actuating_output_indices = actuating_entity_indices
      #     #     else:
      #     #       actuating_output_indices = np.vstack((actuating_output_indices, actuating_entity_indices))
            
      #     #   [[min max][min max][min max]]
        

      #     # actuating_output_indices_list = actuating_output_indices.flatten().tolist()
         
      #     print()
      #     print()
      #     # print('CURRENT ACTUATED CONTROL POINTS SHAPE: ', current_actuated_control_points.shape)
      #     # print('ACTUATING OUTPUT INDICES TYPE: ', type(actuating_output_indices_list))
      #     # print('LENGTH OF ACTUATING OUTPUT INDICES: ', len(actuating_output_indices_list))
      #     # print('ACTUATING OUTPUT INDICES LIST: ', actuating_output_indices_list)
      #     actuated_comp_list = []
      #     for comp_i, component in enumerate(actuation_obj.actuating_components):
      #       min_ind = act_ind_dict[actuation_name][comp_i][0]
      #       max_ind = act_ind_dict[actuation_name][comp_i][1]

      #       actuating_control_points = current_actuated_control_points[min_ind:max_ind+1, :]


      #       # print('CURRENT ACTUATED CONTROL POINTS SHAPE: ', current_actuated_control_points.shape)
      #       # print('Actuating Output Indices: ', type(actuating_output_indices))

      #       # print('Actuating Output Indices flatten shape: ', actuating_output_indices.shape)
      #       # actuating_control_points = self.create_output('actuating_control_points', shape=(num_points,3))

      #       # for i in actuating_output_indices:
      #       #   actuating_control_points = current_actuated_control_points[actuating_output_indices, :]
      #       axis_points = all_points[int(axis.output_starting_ind), :]
      #       origin_points = all_points[int(origin.output_starting_ind), :]


      #       normalize_axis = axis_points / csdl.expand(csdl.pnorm(axis_points), axis_points.shape, 'i->ij')

      #       actuating_points_num = actuating_control_points.shape[:-1]

      #       ''' CURRENTLY IS OKAY BECAUSE WE ACTUATING EVERYTHING RIGHT NOW '''
      #       # actuating_control_points = csdl.reshape(pre_actuated_control_points_expanded[t,:,:], new_shape=(num_points, 3))
      #       print()
      #       print()
      #       print('Actuating points number: ', actuating_points_num)

      #       if actuation_obj.rot_or_trans == 'rot':
      #         quat = self.create_output(f'quat_{act_counter}_{comp_i}_{t}', shape=actuating_points_num + (4,))
              
      #         normalize_axis_x = csdl.reshape(normalize_axis[0,0], new_shape=(1,))
      #         normalize_axis_y = csdl.reshape(normalize_axis[0,1], new_shape=(1,))
      #         normalize_axis_z = csdl.reshape(normalize_axis[0,2], new_shape=(1,))
              
      #         quat[:,0] = csdl.expand(csdl.cos(thetas[t]/2), actuating_points_num + (1,), 'i->ij')
      #         quat[:,1] = csdl.expand(csdl.sin(thetas[t]/2) * normalize_axis_x, actuating_points_num + (1,), 'i->ij')
      #         quat[:,2] = csdl.expand(csdl.sin(thetas[t]/2) * normalize_axis_y, actuating_points_num + (1,), 'i->ij')
      #         quat[:,3] = csdl.expand(csdl.sin(thetas[t]/2) * normalize_axis_z, actuating_points_num + (1,), 'i->ij')
              
      #         origin_points = csdl.reshape(origin_points, new_shape=(3,))
      #         origin_points = csdl.expand(origin_points, actuating_control_points.shape, 'i->ji')

      #         control_points_origin_frame = actuating_control_points - origin_points
              
      #         rotated_control_points_origin_frame = csdl.quatrotvec(quat, control_points_origin_frame)
      #         origin_points = csdl.reshape(origin_points, new_shape=rotated_control_points_origin_frame.shape)
              
      #         joint_actuated_control_points = rotated_control_points_origin_frame + origin_points

      #         # joint_actuated_control_points = csdl.quatrotvec(quat, actuating_control_points)
      #       # elif actuation_obj.rot_or_trans == 'trans':
      #       #   translation_profile = self.declare_variable(f'{actuation_name}_translation_profile_{t}', val=actuation_obj.actuation_profile)

      #       #   # perform actuation
      #       #   joint_actuated_control_points = actuating_control_points + translation_profile

      #       joint_actuated_control_points = csdl.expand(joint_actuated_control_points, (1,) + actuating_points_num + (3,), 'ij->tij')
            
      #       if act_counter == 0:
      #         print(f'HELLO_{comp_i}')
      #         print('COMPONENT NAME:', component.name)
      #         print('MINIMUM IND: ', min_ind)
      #         print('MAXIMUM IND: ', max_ind)
      #         print('ACTUATING CONTROL POINTS SHAPE: ', actuating_control_points.shape)

      #         updated_control_points = self.create_output(f'initial_geometry_{comp_i}_{t}', val=np.reshape(geo.total_cntrl_pts_vector, (1,) + geo.total_cntrl_pts_vector.shape))

      #         updated_control_points[0, min_ind:max_ind+1, :] = joint_actuated_control_points
      #         print('updated control points hash act 0: ', updated_control_points)              
      #         actuated_comp_list.append(updated_control_points)

      #       else:
      #         print('HELLO')
      #         updated_control_points = self.create_output(f'updated_control_points_{act_counter}_{comp_i}_{t}', shape=(1, num_points, 3), val=np.zeros((1, num_points, 3)))
      #         previous_point = outer_act_list[act_counter - 1][comp_i]
      #         # non_actuating_indices = np.delete(range(num_points), actuating_control_points)
      #         if min_ind != 0:
      #           updated_control_points[0, :min_ind ,:] = previous_point[0, :min_ind ,:]
              
      #         if max_ind + 1 < num_points:
      #           updated_control_points[0, max_ind+1: ,:] = previous_point[0, max_ind+1: ,:]

      #         updated_control_points[0, min_ind:max_ind+1, :] = joint_actuated_control_points 
      #         actuated_comp_list.append(updated_control_points)

      #     outer_act_list.append(actuated_comp_list)  
      #     act_counter += 1
            
      #   val_list.append(outer_act_list)
      #   # print('CHECKING TO SEE IF STORES SEPARATE OBJECTS: ', val_list)

      # print('VALUE LIST: ', val_list)
      # actuated_control_points = self.create_output('actuated_control_points', shape=(nt, num_points, 3))
      # for t in range(len(val_list)):
      #   print('t in the val_list loop: ', t)
      #   temp_updated_points = 0
      #   for i in range(len(val_list[t])):
      #     for j in range(len(val_list[t][i])):
      #       print('j in inner for loop: ', j)
      #       temp = val_list[t][i][j]
      #       temp_updated_points = temp + temp_updated_points
        
      #   temp_updated_points = temp_updated_points - ((nt-1) * csdl.reshape(pre_actuated_control_points, new_shape=(1,) + pre_actuated_control_points.shape))
      #   actuated_control_points[t, :, :] = csdl.reshape(temp_updated_points, new_shape=(1,) + pre_actuated_control_points.shape)
      
      # # print('Actuated control points hash: ', actuating_control_points)
      # # for i in range(len(val_list)):
      # #   print('VAL LIST: ', val_list[i])
      # #   actuated_control_points[i, :, :] = val_list[i]


# if __name__ == "__main__":
#     import csdl
#     import numpy as np
#     from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
#     from lsdo_kit.design.design_geometry.tests.csdl_block_update_test import CsdlBlockUpdateTest

#     nxp = 5
#     nyp = 5
#     nzp = 5

#     point000 = np.array([170. ,0. ,100.])
#     point010 = np.array([130., 230., 100.])
#     point001 = np.array([170., 0., 170.])
#     point011 = np.array([130., 230., 170.])
    
#     point100 = np.array([240. ,0. ,100.])
#     point101 = np.array([240. ,0. ,170.])
#     point110 = np.array([200. ,230. ,100.])
#     point111 = np.array([200. ,230. ,170.])

#     control_points = np.zeros((2,2,2,3))
    
#     control_points[0,0,0,:] = point000
#     control_points[0,0,1,:] = point001

#     control_points[0,1,0,:] = point010
#     control_points[0,1,1,:] = point011
    
#     control_points[1,0,0,:] = point100
#     control_points[1,0,1,:] = point101
    
#     control_points[1,1,0,:] = point110
#     control_points[1,1,1,:] = point111

#     ffd_control_points = create_ffd(control_points, nxp, nyp, nzp)

#     # Create ffd_blocks list, and the csdl variables needed

#     ffd_blocks = [ffd_control_points]

#     sim = Simulator(CsdlBlockUpdateTest())
#     sim.run()
#     sim.visualize_model()
