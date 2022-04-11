import numpy as np

from openmdao.api import ExplicitComponent
'''
Quaternions: real -> imaginary (left-right) or (top-down)
'''


class ComputeQuatRateComp(ExplicitComponent):
    def initialize(self):
        self.options.declare('angvelshape', types=tuple)
        self.options.declare('quatshape', types=tuple)
        self.options.declare('axis', types=int)
        self.options.declare('angvelval', types=np.ndarray)
        self.options.declare('quatval', types=np.ndarray)

        self.options.declare(
            'angvel_name',
            types=str,
            desc="This the name of rotation axis for the quaternion")

        self.options.declare(
            'quat_name',
            types=str,
            desc="This the name of rotation axis for the quaternion")

        self.options.declare('out_name',
                             types=str,
                             desc="The name of the quaternion")

    def setup(self):
        axis = self.options['axis']

        angvelshape = self.options['angvelshape']
        angvelval = self.options['angvelval']
        angvel_name = self.options['angvel_name']

        quatshape = self.options['quatshape']
        quatval = self.options['quatval']
        quat_name = self.options['quat_name']

        out_name = self.options['out_name']

        # Defining the shape of the cross matrix needed for the multiplication
        # shape_list = list(angvelshape)
        # rate_shape = tuple(np.delete(shape_list, axis)) + (4, 4)
        # # print(rate_shape)
        # self.ratemat = np.zeros(rate_shape)

        # Adding the angular velocity as an input
        self.add_input(angvel_name, shape=angvelshape, val=angvelval)

        # Defining the quaternion vector for multiplication
        shape_list = list(quatshape)
        quat_shape = (tuple(np.delete(shape_list, axis))) + (4, 1)
        print(quat_shape)
        self.quatvec = np.zeros(quat_shape)

        # Adding the quaternion input
        self.add_input(quat_name, shape=quatshape, val=quatval)

        self.add_output(out_name, shape=quat_shape)

        # -------------- PARTIALS -------------------

        # rows = np.tile(np.arange(np.prod(shape)), 4)
        # cols = np.repeat(np.arange(4), 4)

        # self.declare_partials(out_name, in_name, rows=, cols=cols)
        self.declare_partials(of='*', wrt='*', method='cs')

    def compute(self, inputs, outputs):
        axis = self.options['axis']
        angvel_name = self.options['angvel_name']
        angvelshape = self.options['angvelshape']
        quatshape = self.options['quatshape']
        quat_name = self.options['quat_name']
        out_name = self.options['out_name']

        # --------- SETTING UP CROSS MATRIX BASED ON ANGULAR VELOCITY ----------
        ang_x = np.real(np.take(inputs[angvel_name], 0, axis=axis))
        ang_y = np.real(np.take(inputs[angvel_name], 1, axis=axis))
        ang_z = np.real(np.take(inputs[angvel_name], 2, axis=axis))

        quat_r = np.real(np.take(inputs[quat_name], 0, axis=axis))
        quat_i = np.real(np.take(inputs[quat_name], 1, axis=axis))
        quat_j = np.real(np.take(inputs[quat_name], 2, axis=axis))
        quat_k = np.real(np.take(inputs[quat_name], 3, axis=axis))

        shape_list = list(angvelshape)
        rate_shape = tuple(np.delete(shape_list, axis)) + (4, 4)
        ratemat = np.zeros(rate_shape)

        shape_list = list(quatshape)
        quat_shape = (tuple(np.delete(shape_list, axis))) + (4, 1)
        quatvec = np.zeros(quat_shape)

        # outputs[out_name] = np.put_along_axis(outputs[out_name], )

        ratemat[..., 0, 1] = -ang_x
        ratemat[..., 0, 2] = -ang_y
        ratemat[..., 0, 3] = -ang_z
        ratemat[..., 1, 0] = ang_x
        ratemat[..., 2, 0] = ang_y
        ratemat[..., 3, 0] = ang_z

        ratemat[..., 1, 2] = -ang_z
        ratemat[..., 1, 3] = ang_y
        ratemat[..., 2, 1] = ang_z
        ratemat[..., 2, 3] = -ang_x
        ratemat[..., 3, 1] = -ang_y
        ratemat[..., 3, 2] = ang_x

        # print(ratemat.shape)
        # -------- Setting up Quat vector based on quaternion -------
        # print(quatvec.shape)

        quatvec[..., 0, 0] = quat_r
        quatvec[..., 1, 0] = quat_i
        quatvec[..., 2, 0] = quat_j
        quatvec[..., 3, 0] = quat_k

        outputs[out_name] = np.matmul(ratemat, quatvec)


if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    angvelshape = (3, 5)

    quatshape = (4, 5)

    axis = 0

    prob = Problem()

    angvelval = np.prod(angvelshape)
    quatval = np.prod(quatshape)

    comp = IndepVarComp()

    comp.add_output('angvel', val=np.arange(angvelval).reshape(angvelshape))
    comp.add_output('quat', val=np.arange(quatval).reshape(quatshape))

    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = ComputeQuatRateComp(
        angvelshape=angvelshape,
        quatshape=quatshape,
        angvelval=np.arange(angvelval).reshape(angvelshape),
        quatval=np.arange(quatval).reshape(quatshape),
        axis=axis,
        angvel_name='angular_vel',
        quat_name='quat',
        out_name='compute_quat_rate',
    )
    prob.model.add_subsystem('compute_quat_rate_comp', comp, promotes=['*'])

    prob.setup()
    prob.run_model()
    # prob.check_partials(compact_print=True)

    # prob.model.list_inputs(print_arrays=True)
    # prob.model.list_outputs(print_arrays=True)