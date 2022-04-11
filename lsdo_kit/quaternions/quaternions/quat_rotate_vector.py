import numpy as np
# import csdl
from openmdao.api import ExplicitComponent
import os
os.chdir("../lsdo_geo/quaternions/quaternions")

'''
Quaternions: real -> imaginary (left-right) or (top-down)
'''


class QuatRotateVectorComp(ExplicitComponent):
    def initialize(self):
        # self.options.declare('quatshape', types=tuple)
        # self.options.declare('vecshape', types=tuple)
        self.options.declare('quat_vals', types=np.ndarray)
        self.options.declare('vec_vals', types=np.ndarray)

        self.options.declare('shape', types=tuple) 

        self.options.declare(
            'quat_name',
            types=str,
            desc="This the name of rotation axis for the quaternion")

        self.options.declare(
            'vec_name',
            types=str,
            desc="This the name of angle which will be rotated about the axis")

        self.options.declare('out_name',
                             types=str,
                             desc="The name of the quaternion")

    def setup(self):

        quat_vals = self.options['quat_vals']
        vec_vals = self.options['vec_vals']
        shape = self.options['shape']

        quat_name = self.options['quat_name']
        vec_name = self.options['vec_name']
        out_name = self.options['out_name']

        self.add_input(quat_name, shape=shape + (4,), val=quat_vals)
        self.add_input(vec_name, shape=shape + (3,), val=vec_vals)
        self.add_output(out_name, shape=shape + (3,))

        ''' Formatting shapes of output vectors '''
        size = np.prod(shape)
        out_indices = np.arange(size*3).reshape(shape +(3,))
        vec_indices = np.arange(size*3).reshape(shape +(3,))
        quat_indices = np.arange(size*4).reshape(shape +(4,))

        r0 = np.zeros(shape + (3,3))
        c0 = np.zeros(shape + (3,3))
        
        ''' Creating rows and columns for partials'''
        for i in range(3):
            for j in range(3):
                r0[...,i,j] = out_indices[...,i]
                c0[...,i,j] = vec_indices[...,j]
                
        self.declare_partials(out_name, vec_name, rows=r0.flatten(), cols=c0.flatten())
        
        r1 = np.zeros(shape + (3,4))
        c1 = np.zeros(shape + (3,4))

        for i in range(3):
            for j in range(4):
                r1[...,i,j] = out_indices[...,i]
                c1[...,i,j] = quat_indices[...,j]

        self.declare_partials(out_name, quat_name, rows=r1.flatten(), cols=c1.flatten())

    def dot_quat_vec(self, quat, vec):
        dot_quat_vec = vec[...,0]*quat[...,1] + vec[...,1]*quat[...,2] + vec[...,2]*quat[...,3]
        return dot_quat_vec

    def dot_quat_quat(self, quat):
        dot_quat_quat = quat[...,1]**2 + quat[...,2]**2 + quat[...,3]**2
        return dot_quat_quat

    def execute(self, quat, vec):
        # print('------------------')
        # print('Quaternion xyz: ', quat[...,1:4])
        # print('Quaternion xyz shape: ', quat[...,1:4].shape)
        # print('Quaternion real: ', quat[...,0])
        # print('Quaternion real shape: ', quat[...,0].shape)

        # temp1 = np.einsum('ij,ij->i', quat[...,1:4], vec)
        # temp1dot = np.dot(quat[0,1:4], vec[0,:])
        # temp1dot2 = np.dot(quat[1,1:4], vec[1,:])
        
        # print('First Dot Product: ', temp1)

        # temp2 = (quat[...,1:4].T * temp1).T
        # temp2test = quat[0,1:4] * temp1[0]
        # temp2test1 = quat[1,1:4] * temp1[1]
        # print('Multiplied First Dot Product by QuatVector: ', temp2)

        # temp3 = np.einsum('ij,ij->i', quat[...,1:4], quat[...,1:4])
        # temp3test = np.dot(quat[0,1:4], quat[0,1:4])
        # temp3test1 = np.dot(quat[1,1:4], quat[1,1:4])
        # print('Dot product between quat and quat imag comp: ', temp3)

        # temp4 = (quat[...,0] ** 2) - temp3
        # temp4test = quat[0,0]**2 - temp3[0]
        # temp4test1 = quat[1,0]**2 - temp3[1]
        # print('Quat Real Squared - DotProd: ', temp4)

        # temp5 = (vec.T * temp4).T
        # temp5test = temp4[0] * vec[0,:]
        # temp5test1 = temp4[1] * vec[1,:]
        # print('Quat Real Squared - DotProd multiplied by Vector: ', temp5)

        # temp6 = np.cross(quat[...,1:4], vec)
        # temp6test = np.cross(quat[0,1:4], vec[0,:])
        # temp6test1 = np.cross(quat[1,1:4], vec[1,:])
        # print('Cross Product between the Quaternion and Vector: ', temp6)

        # temp7 = 2*(quatval[...,0] * temp6.T).T
        # temp7test = 2*quat[0,0] * np.cross(quat[0,1:4], vec[0,:])
        # temp7test1 = 2*quat[1,0] * np.cross(quat[1,1:4], vec[1,:])
        # print('Cross Product Multiplied by Quaternion Real, Multiplied by 2: ', temp7)

        # temp8=2.0 * (quat[...,1:4].T * np.einsum('...ij,...ij->...i', quat[...,1:4], vec)).T \
        #     + (vec.T * ((quat[...,0] ** 2) - np.einsum('...ij,...ij->...i', quat[...,1:4], quat[...,1:4]))).T \
        #     + 2.0 * (np.cross(quat[...,1:4], vec).T * quatval[...,0]).T
        # print('Final Result: ', temp8)


        # temp1 = (quat[...,1:4].T * np.einsum('ij,ij->i', quat[...,1:4], vec)).T
        # temp2 = (vec.T * ((quat[...,0] ** 2) - np.einsum('ij,ij->i', quat[...,1:4], quat[...,1:4]))).T      

        # print('Temp1: ', temp1)
        # print('Temp2: ', temp2)
        # print('cross products: ', (np.cross(quat[...,1:4], vec).T * quatval[...,0]).T)

        # print('quat1-4: ', quat[...,1:4].shape)
        shape = self.options['shape']   
        temp = np.zeros(shape + (3,))
        temp[...,0] = 2 * self.dot_quat_vec(quat,vec) * quat[...,1] + (quat[...,0]**2 - self.dot_quat_quat(quat)) * vec[...,0] + 2*quat[...,0] * (quat[...,2]*vec[...,2] - quat[...,3]*vec[...,1])  
        
        temp[...,1] = 2 * self.dot_quat_vec(quat,vec) * quat[...,2] + (quat[...,0]**2 - self.dot_quat_quat(quat)) * vec[...,1] + 2*quat[...,0] * (quat[...,3]*vec[...,0] - quat[...,1]*vec[...,2])

        temp[...,2] = 2 * self.dot_quat_vec(quat,vec) * quat[...,3] + (quat[...,0]**2 - self.dot_quat_quat(quat)) * vec[...,2] + 2*quat[...,0] * (quat[...,1]*vec[...,1] - quat[...,2]*vec[...,0])


        # test = 2 * quat[...,1:4] * np.expand_dims(np.einsum('...i,...i->...', quat[...,1:4], vec), axis=-1)

        # print('quat0: ', quat[...,0].shape)

        # test2 = vec * np.expand_dims((quat[...,0] ** 2) - np.einsum('...i,...i->...', quat[...,1:4], quat[...,1:4]),axis=-1)

        # test3 = 2 * np.expand_dims(quatval[...,0],axis=-1) * np.cross(quat[...,1:4], vec)

        # out = test + test2 + test3

        # out = 2.0 * (quat[...,1:4].T * np.einsum('...i,...i->...', quat[...,1:4], vec)).T \
        #     + (vec.T * ((quat[...,0] ** 2) - np.einsum('...i,...i->...', quat[...,1:4], quat[...,1:4]))).T \
        #     + 2.0 * (quatval[...,0] * np.cross(quat[...,1:4], vec).T).T

        # print('Output: ', out)
        # print('Output: ', out.shape)
        
        return temp

    def compute(self, inputs, outputs):
        quat_name = self.options['quat_name']
        vec_name = self.options['vec_name']
        out_name = self.options['out_name']

        outputs[out_name] = self.execute(inputs[quat_name], inputs[vec_name])

    def compute_partials(self, inputs, partials):
        quat_name = self.options['quat_name']
        vec_name = self.options['vec_name']
        out_name = self.options['out_name']


        quat = inputs[quat_name]
        vec = inputs[vec_name]

        shape = self.options['shape']

        temp = np.zeros(shape + (3,3))

        temp[...,0,0] = 2*quat[...,1]**2 + quat[...,0]**2 - quat[...,1]**2 - quat[...,2]**2 - quat[...,3]**2 
        temp[...,0,1] = 2*quat[...,2]*quat[...,1] - 2*quat[...,0]*quat[...,3]
        temp[...,0,2] = 2*quat[...,3]*quat[...,1] + 2*quat[...,0]*quat[...,2]
        
        temp[...,1,0] = 2*quat[...,1]*quat[...,2] + 2*quat[...,0]*quat[...,3]
        temp[...,1,1] = 2*quat[...,2]**2 + quat[...,0]**2 - quat[...,1]**2 - quat[...,2]**2 - quat[...,3]**2 
        temp[...,1,2] = 2*quat[...,3]*quat[...,2] - 2*quat[...,0]*quat[...,1]

        temp[...,2,0] = 2*quat[...,1]*quat[...,3] - 2*quat[...,0]*quat[...,2]
        temp[...,2,1] = 2*quat[...,2]*quat[...,3] + 2*quat[...,0]*quat[...,1]
        temp[...,2,2] = 2*quat[...,3]**2 + quat[...,0]**2 - quat[...,1]**2 - quat[...,2]**2 - quat[...,3]**2

        partials[out_name, vec_name] = temp.flatten()

        temp1 = np.zeros(shape + (3,4))

        temp1[...,0,0] = 2*quat[...,0]*vec[...,0] + 2*quat[...,2]*vec[...,2] - 2*quat[...,3]*vec[...,1]
        temp1[...,0,1] = 4*vec[...,0]*quat[...,1] + 2*vec[...,1]*quat[...,2] + 2*vec[...,2]*quat[...,3] - 2*vec[...,0]*quat[...,1] 
        temp1[...,0,2] = 2*vec[...,1]*quat[...,1] - 2*quat[...,2]*vec[...,0] + 2*quat[...,0]*vec[...,2]
        temp1[...,0,3] = 2*vec[...,2]*quat[...,1] - 2*quat[...,3]*vec[...,0] - 2*quat[...,0]*vec[...,1]

        temp1[...,1,0] = 2*quat[...,0]*vec[...,1] + 2*quat[...,3]*vec[...,0] - 2*quat[...,1]*vec[...,2]
        temp1[...,1,1] = 2*vec[...,0]*quat[...,2] - 2*quat[...,1]*vec[...,1] - 2*quat[...,0]*vec[...,2]
        temp1[...,1,2] = 4*vec[...,1]*quat[...,2] + 2*vec[...,0]*quat[...,1] + 2*vec[...,2]*quat[...,3] - 2*quat[...,2]*vec[...,1]  
        temp1[...,1,3] = 2*vec[...,2]*quat[...,2] - 2*quat[...,3]*vec[...,1] + 2*quat[...,0]*vec[...,0]

        temp1[...,2,0] = 2*quat[...,0]*vec[...,2] + 2*quat[...,1]*vec[...,1] - 2*quat[...,2]*vec[...,0]
        temp1[...,2,1] = 2*vec[...,0]*quat[...,3] - 2*quat[...,1]*vec[...,2] + 2*quat[...,0]*vec[...,1]
        temp1[...,2,2] = 2*vec[...,1]*quat[...,3] - 2*quat[...,2]*vec[...,2] - 2*quat[...,0]*vec[...,0]
        temp1[...,2,3] = 2*vec[...,0]*quat[...,1] + 2*vec[...,1]*quat[...,2] + 4*vec[...,2]*quat[...,3] - 2*quat[...,3]*vec[...,2] 
         
        partials[out_name, quat_name] = temp1.flatten()
    
if __name__ == "__main__":
    from openmdao.api import Problem, IndepVarComp

    quatshape = (2,4)
    vecshape = (2,3)
    prob = Problem()

    # quatval = np.array([[0.18257419, 0.36514837, 0.54772256, 0.73029674], [0.04724383, 0.63094826, 0.42525845, 0.64716889]])
    quatval = np.ones(quatshape)
    quatval = np.reshape(np.arange(np.prod(quatshape)), (quatshape))
    vecval = np.reshape(np.arange(np.prod(vecshape)), vecshape)

    comp = IndepVarComp()
    comp.add_output('quat', val=quatval)
    comp.add_output('vec', val=vecval)

    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = QuatRotateVectorComp(
        shape=(2,),
        quat_vals=quatval,
        vec_vals=vecval,
        quat_name='quat',
        vec_name='vec',
        out_name='output',
    )
    prob.model.add_subsystem('rotated_vec', comp, promotes=['*'])

    prob.setup()
    prob.run_model()
    prob.check_partials(compact_print=True)
    prob.model.list_inputs(print_arrays=True)
    prob.model.list_outputs(print_arrays=True)













 
    # # prob = Problem()

    # quatval = np.random.random(quatshape)
    # vecval = np.random.random(vecshape)

    # prob.setup()
    # # prob.run_model()


    # comp = IndepVarComp()
    # comp.add_output('quat', val=quatval)
    # comp.add_output('vec', val=vecval)

    # prob.model.add_subsystem('ivc', comp, promotes=['*'])

    # comp = QuatRotateVectorComp(
    #     quatshape=quatshape,
    #     vecshape=vecshape,
    #     quatval=quatval,
    #     vecval=vecval,
    #     axis=axis,
    #     quat_name='quat',
    #     vec_name='vec',
    #     out_name='rotvec',
    # )
    # prob.model.add_subsystem('rotated_vec', comp, promotes=['*'])


    # prob.model.list_inputs(print_arrays=True)
    # prob.model.list_outputs(print_arrays=True)

    #     def get_indices(shape):
    #         return np.arange(np.prod(shape)).reshape(shape)

    #     self.vectorized_shape = list(quatshape)
    #     self.vectorized_shape.pop(axis)

    #     quat_indices = get_indices(quat_shape)
    #     vec_indices = get_indices(vec_shape)
    #     out_indices = get_indices(out_shape)

    #     rows = np.zeros(shape + (3, 4), int)
    #     for ind in range(4):
    #         rows[..., :, ind] = out_indices
    #     cols = np.zeros(shape + (3, 4), int)
    #     for ind in range(3):
    #         cols[..., ind, :] = quat_indices
    #     self.declare_partials(out_name, quat_name, rows=rows, cols=cols, method='cs')
            
    #     rows = np.zeros(shape + (3, 3), int)
    #     for ind in range(3):
    #         rows[..., :, ind] = out_indices
    #     cols = np.zeros(shape + (3, 3), int)
    #     for ind in range(3):
    #         cols[..., ind, :] = vec_indices
    #     self.declare_partials(out_name, vec_name, rows=rows, cols=cols, method='cs')

    # def execute(self, quat, vec):
    #     axis = self.options['axis']
    #     s = np.expand_dims(np.take(quat, 0, axis=axis), axis=axis)

    #     # index = [1,2,3] because we assume that the i,j,k components follow the 0th index
    #     qxyz = np.take(quat, [1, 2, 3], axis=axis)

    #     dotproduv = np.sum(qxyz * vec, axis=axis, keepdims=True)
    #     dotproduu = np.sum(qxyz * qxyz, axis=axis, keepdims=True)

    #     out = 2.0 * dotproduv * qxyz + (
    #         s * s - dotproduu) * vec + 2.0 * s * np.cross(qxyz, vec, axis=axis)

    #     return out

    # def compute(self, inputs, outputs):
    #     quat_name = self.options['quat_name']
    #     vec_name = self.options['vec_name']
    #     out_name = self.options['out_name']

    #     outputs[out_name] = self.execute(inputs[quat_name], inputs[vec_name])
       
    # def compute_partials(self, inputs, partials):
    #     quat_name = self.options['quat_name']
    #     vec_name = self.options['vec_name']
    #     out_name = self.options['out_name']

    #     h = 1e-20
    #     ih = complex(0,h)

    #     quat = np.array(inputs[quat_name])
    #     vec  = np.array(inputs[vec_name])

    #     partial_quat = partials[out_name, quat_name].reshape(self.vectorized_shape + (3, 4))
    #     partial_vec  = partials[out_name,  vec_name].reshape(self.vectorized_shape + (3, 3))

    #     for ind in range(4):
    #         quat[..., 0] += ih
    #         out = self.execute(quat, vec)
    #         quat[..., 0] -= ih
    #         partial_quat[:, :, :, :, ind] = out.imag / h

    #     for ind in range(3):
    #         vec[..., 0] += ih
    #         out = self.execute(quat, vec)
    #         vec[..., 0] -= ih
    #         partial_vec[:, :, :, :, ind] = out.imag / h

