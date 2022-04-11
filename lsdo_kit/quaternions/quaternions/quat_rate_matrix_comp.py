import numpy as np

from openmdao.api import ExplicitComponent
'''
Quaternions: real -> imaginary (left-right) or (top-down)
'''


class QuatRateMatrixComp(ExplicitComponent):
    def initialize(self):
        self.options.declare('angvelshape', types=tuple)
        self.options.declare('axis', types=int)
        self.options.declare('angvelval', types=np.ndarray)

        self.options.declare(
            'angvel_name',
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
        out_name = self.options['out_name']

        self.add_input(angvel_name, shape=angvelshape, val=angvelval)

        shape_list = list(angvelshape)
        self.out_shape = (np.prod(np.delete(shape_list, axis)), 4, 4)

        self.ratemat = np.zeros(self.out_shape)

        self.add_output(out_name, shape=self.out_shape)

        # rows = np.tile(np.arange(np.prod(shape)), 4)
        # cols = np.repeat(np.arange(4), 4)

        # self.declare_partials(out_name, in_name, rows=, cols=cols)
        self.declare_partials(of='*', wrt='*', method='cs')

    def compute(self, inputs, outputs):
        axis = self.options['axis']
        angvel_name = self.options['angvel_name']
        out_name = self.options['out_name']

        vx = np.take(inputs[angvel_name], 0, axis=axis).flatten()
        vy = np.take(inputs[angvel_name], 1, axis=axis).flatten()
        vz = np.take(inputs[angvel_name], 2, axis=axis).flatten()

        outputs[out_name][..., 0, 1] = -vx
        outputs[out_name][..., 0, 2] = -vy
        outputs[out_name][..., 0, 3] = -vz
        outputs[out_name][..., 1, 0] = vx
        outputs[out_name][..., 2, 0] = vy
        outputs[out_name][..., 3, 0] = vz

        outputs[out_name][..., 1, 2] = -vz
        outputs[out_name][..., 1, 3] = vy
        outputs[out_name][..., 2, 1] = vz
        outputs[out_name][..., 2, 3] = -vx
        outputs[out_name][..., 3, 1] = -vy
        outputs[out_name][..., 3, 2] = vx


if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    angvelshape = (2, 4, 3, 5)
    axis = 2

    prob = Problem()

    angvelval = np.prod(angvelshape)

    comp = IndepVarComp()
    comp.add_output('angvel', val=np.arange(angvelval).reshape(angvelshape))

    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = QuatRateMatrixComp(
        angvelshape=angvelshape,
        angvelval=np.arange(angvelval).reshape(angvelshape),
        axis=axis,
        angvel_name='angvel',
        out_name='quatratemat',
    )
    prob.model.add_subsystem('quatratemat_comp', comp, promotes=['*'])

    prob.setup()
    prob.run_model()
    # prob.check_partials(compact_print=True)

    # prob.model.list_inputs(print_arrays=True)
    # prob.model.list_outputs(print_arrays=True)