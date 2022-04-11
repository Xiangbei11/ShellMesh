import csdl

class EvaluatePointsets(csdl.Model):
    '''
        This will not run if no pointsets exist in DesignGeometry

        REMINDER:
            - Check offset_map 
            - Discuss what to do when no point_sets exist
    '''
    def initialize(self):
        self.parameters.declare('design_geometry_obj')

    def define(self):
        design_geometry_obj = self.parameters['design_geometry_obj']
        design_geometry_obj.assemble()
        eval_map = self.create_input('eval_map', val=design_geometry_obj.eval_map.todense())
        offset_map = self.create_input('offset_map', val=design_geometry_obj.offset_eval_map)
        control_points = self.declare_variable('control_points', val=design_geometry_obj.total_cntrl_pts_vector)

        # print('Offset_map: ', design_geometry_obj.offset_eval_map.shape)

        points = csdl.matmat(eval_map, control_points)

        # print('Control Points Vector Shape: ', design_geometry_obj.total_cntrl_pts_vector.shape)
        # print('Eval Map Shape: ', design_geometry_obj.eval_map.shape)
        
        self.register_output('points', points)
        # updated_points = control_points + offset_map        

# def evaluate_pointsets(design_geometry_obj, application_control_points):
#     design_geometry_obj.assemble()
#     eval_map = design_geometry_obj.eval_map.todense()
#     offset_map = design_geometry_obj.offset_eval_map

#     points = csdl.matmat(eval_map, application_control_points)
#     updated_points = application_control_points + offset_map
    
#     return updated_points