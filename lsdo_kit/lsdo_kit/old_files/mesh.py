import numpy as np

class Mesh:
    '''
    Class for representing meshes to be fed into analysis components.
    '''

    def __init__(self, name, pointset_list=[]) -> None:
        self.name = name
        self.pointset_list = []
        self.features_list = []
        self.num_mesh_points = 0

        for pointset in pointset_list:
            self.add_pointset(pointset)

    def add_pointset(self, pointset, name=None):
        '''
        Create/add on to the single instance of evaluated_pointsets
        '''
        if name is None:
            name = f'pointset_{pointset.id}'
        pointset.name = name

        self.pointset_list.append(pointset)
        self.num_mesh_points += np.cumprod(pointset.shape)[-2]

    def add_feature(self, feature):

        self.features_list.append(feature)



# class BSplineMesh(Mesh):
#     '''
#     Class for representing bspline meshes to be fed into analysis components like IMGA.
#     '''

#     def __init__(self) -> None:
#         pass

#     def add_pointset(self, pointset, name=None):
#         '''
#         Create/add on to the single instance of evaluated_pointsets
#         '''
#         if name is None:
#             name = f'pointset_{pointset.id}'

#         self.pointset_list[name] = pointset

#     def fit_bspline_entities(self, points):
#         '''Least square fit b-spline surface'''
#         entity_starting_point = 0
#         i = 0
#         for pointer_id in self.evaluated_pointsets.parent_pointers:
#             component_shape = self.pointsets_dict[pointer_id].shape
#             #print(component_shape)
#             #print(component_shape[:-1])
#             #print(entity_starting_point)
#             if len(component_shape[:-1]) == 0:  # is point
#                 print('fitting points has not been implemented yet')
#                 pass        #is point
#             elif len(component_shape[:-1]) == 1:  # is curve
#                 if i==0:
#                     entity_points_be_fitted = points[0:component_shape[0],:]
#                 else:
#                     entity_points_be_fitted = points[entity_starting_point:(entity_starting_point+component_shape[0]),:]
#                 entity_starting_point += component_shape[0]
#                 curve = fitting.approximate_curve(entity_points_be_fitted, 3)   #TODO hardcoded cubic bspline
#                 bspline_entity_curve = BSplineCurve(
#                     name=self.registered_names[i],
#                     order_u=curve.order,
#                     shape=np.array[len(curve.ctrlpts), 3],
#                     control_points=curve.ctrlpts,
#                     knots_u=curve.knotvector)
#                 self.output_bspline_entity_dict[self.registered_names[i]] = bspline_entity_curve
#             elif len(component_shape[:-1]) == 2:  # is surface
#                 if i==0:
#                     num_pts = component_shape[0] * component_shape[1]
#                     entity_points_be_fitted = points[0:num_pts,:]
#                 else:
#                     num_pts = component_shape[0] * component_shape[1]
#                     entity_points_be_fitted = points[entity_starting_point:(entity_starting_point+num_pts),:]
#                 entity_starting_point += num_pts
#                 order_u_fitted = 4
#                 order_v_fitted = 4
#                 num_control_points_u_fitted = component_shape[0] - 1
#                 num_control_points_v_fitted = component_shape[1] - 1
#                 num_points_u_fitted = component_shape[0]
#                 num_points_v_fitted = component_shape[1]

#                 nnz = num_points_u_fitted * num_points_v_fitted * order_u_fitted * order_v_fitted
#                 data = np.zeros(nnz)
#                 row_indices = np.zeros(nnz, np.int32)
#                 col_indices = np.zeros(nnz, np.int32)
#                 u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u_fitted), np.ones(num_points_v_fitted)).flatten()
#                 v_vec = np.einsum('i,j->ij', np.ones(num_points_u_fitted), np.linspace(0., 1., num_points_v_fitted)).flatten()
#                 knot_vector_u = np.zeros(num_control_points_u_fitted+order_u_fitted)
#                 knot_vector_v = np.zeros(num_control_points_v_fitted+order_v_fitted)
#                 get_open_uniform(order_u_fitted, num_control_points_u_fitted, knot_vector_u)
#                 get_open_uniform(order_v_fitted, num_control_points_v_fitted, knot_vector_v)
#                 get_basis_surface_matrix(
#                     order_u_fitted, num_control_points_u_fitted, 0, u_vec, knot_vector_u,
#                     order_v_fitted, num_control_points_v_fitted, 0, v_vec, knot_vector_v,
#                     num_points_u_fitted * num_points_v_fitted, data, row_indices, col_indices,
#                 )
#                 basis0 = sps.csc_matrix(
#                     (data, (row_indices, col_indices)), 
#                     shape=(num_points_u_fitted * num_points_v_fitted, num_control_points_u_fitted * num_control_points_v_fitted),
#                 )

#                 a = np.matmul(basis0.toarray().T, basis0.toarray())
#                 if np.linalg.det(a) == 0:
#                     cps_fitted,_,_,_ = np.linalg.lstsq(a, np.matmul(basis0.toarray().T, entity_points_be_fitted), rcond=None)
#                 else: 
#                     cps_fitted = np.linalg.solve(a, np.matmul(basis0.toarray().T, entity_points_be_fitted))            
            
#                 bspline_entity_surface = BSplineSurface(
#                     name=self.registered_names[i],
#                     order_u=order_u_fitted,
#                     order_v=order_v_fitted,
#                     shape=np.array([num_control_points_u_fitted, num_control_points_v_fitted, 3]),
#                     control_points=np.array(cps_fitted).reshape((num_control_points_u_fitted*num_control_points_v_fitted,3)),
#                     knots_u=np.array(knot_vector_u),
#                     knots_v=np.array(knot_vector_v))
#                 self.output_bspline_entity_dict[self.registered_names[i]] = bspline_entity_surface

#             elif len(component_shape[:-1]) == 3:  # is volume
#                 print('fitting BSplineVolume has not been implemented yet')
#                 pass
#             i += 1
#         pass 