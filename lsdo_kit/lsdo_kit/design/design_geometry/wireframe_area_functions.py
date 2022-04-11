import numpy as np

def starting_indices(n,m):
    '''
    Objective: When given a mesh of points, it gives the upper left coordinates of each
    rectangular mesh.

    Input: Dimensions of mesh

    Output: Upper left coordinate of a rectangular mesh
    '''
    index_list = []
    for i in range(n-1):
        for j in range(m-1):
            temp=[i,j]
            index_list.append(temp)
    return index_list

def get_area(mesh_coordinates):
  '''
  Objective: Get area of a mesh quadrilateral from a cross product operation

  Input: Coordinates of each mesh quadrilateral in a list

  Output: Area of mesh
  '''
  up_vectors = []
  side_vectors = []
  area = 0
  for i in range(len(mesh_coordinates)):
    up_vectors.append(mesh_coordinates[i][0]-mesh_coordinates[i][2])
    side_vectors.append(mesh_coordinates[i][2]-mesh_coordinates[i][3])
  cross_list = np.cross(up_vectors,side_vectors)
  for i in range(len(cross_list)):
    area += np.linalg.norm(cross_list[i])
  return area
  
def get_vectors(mesh_coordinates):
    '''
    Objective: Get area of a mesh quadrilateral from a cross product operation

    Input: Coordinates of each mesh quadrilateral in a list

    Output: Area of mesh
    '''
    up_vectors = []
    side_vectors = []
    area = 0
    for i in range(len(mesh_coordinates)):
        up_vectors.append(mesh_coordinates[i][0]-mesh_coordinates[i][2])
        side_vectors.append(mesh_coordinates[i][2]-mesh_coordinates[i][3])
    return up_vectors, side_vectors

def mesh_coordinates_list(n,m,matrix):
  ''' 
  Objective: Collects the vertices of each mesh rectangle into an array
  
  Inputs:

  n = number of rows in matrix
  m = number of columns in matrix
  matrix = matrix of coordinate points, essentially a point cloud

  Outputs:

  Array of each mesh rectangle's vertex points.
  
  '''

  original_matrix_shape = (int(np.sqrt(matrix.shape[0])),int(np.sqrt(matrix.shape[0])))
  mesh_coordinates = []
  indexes = starting_indices(n,m)
  for i in range(len(starting_indices(n,m))):
    temp = []
    for j in range(0,4):
      start_x_index = indexes[i][0]
      start_y_index = indexes[i][1] 
      if j == 0:
        first_index = matrix[np.ravel_multi_index((start_x_index,start_y_index),original_matrix_shape)]
        temp.append(first_index)
      elif j == 1:
        start_y_index += 1
        second_index = matrix[np.ravel_multi_index((start_x_index,start_y_index),original_matrix_shape)]
        temp.append(second_index)
      elif j == 2:
        start_x_index += 1
        third_index = matrix[np.ravel_multi_index((start_x_index,start_y_index),original_matrix_shape)]
        temp.append(third_index)
      elif j == 3:
        start_x_index += 1
        start_y_index += 1
        fourth_index = matrix[np.ravel_multi_index((start_x_index,start_y_index),original_matrix_shape)]
        temp.append(fourth_index)
    mesh_coordinates.append(temp)
  return mesh_coordinates