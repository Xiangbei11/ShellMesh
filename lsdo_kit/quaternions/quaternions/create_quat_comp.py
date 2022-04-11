import numpy as np

from openmdao.api import ExplicitComponent

from csdl.utils.get_array_indices import get_array_indices

alphabet = 'abcdefghijklmnopqrstuvw'
'''
Quaternions: real -> imaginary (left-right) or (top-down)
'''


class CreateQuatComp(ExplicitComponent):
    def initialize(self):
        self.options.declare('vecshape', types=tuple)
        self.options.declare('angleshape', types=tuple)
        self.options.declare('axis', types=int)
        self.options.declare('vecval', types=np.ndarray)
        self.options.declare('angleval', types=np.ndarray)

        self.options.declare(
            'vec_name',
            types=str,
            desc="This the name of rotation axis for the quaternion")

        self.options.declare(
            'angle_name',
            types=str,
            desc="This the name of angle which will be rotated about the axis")

        self.options.declare('out_name',
                             types=str,
                             desc="The name of the quaternion")

    def setup(self):
        vecshape = self.options['vecshape']
        angleshape = self.options['angleshape']
        axis = self.options['axis']
        vecval = self.options['vecval']
        angleval = self.options['angleval']
        vec_name = self.options['vec_name']
        angle_name = self.options['angle_name']
        out_name = self.options['out_name']

        self.add_input(vec_name, shape=vecshape, val=vecval)

        # ANGLE DIMENSION must have a 1 in the specified axis. vec shape = (2,3,4), angle shape = (2,1,4)
        self.add_input(angle_name, shape=angleshape, val=angleval)

        ''' Defining what the output shape should be '''
        shape_list = list(vecshape)
        shape_list[axis] = 4
        out_shape = tuple(shape_list)

        self.add_output(out_name, shape=out_shape)

        shape_list = list(vecshape)
        shape_list[axis] = 1
        quat_real_shape = tuple(shape_list)

        self.quat_real = np.ones(quat_real_shape)

        # Partials rows and cols wrt angle
        indices = get_array_indices(*angleshape)

        # self.shape_without_axis = shape[:axis] + shape[axis + 1:]

        # ones = np.ones(1, int)
        # rank = len(self.shape_without_axis)

        # einsum_string_rows = '{}y{},z->{}{}yz'.format(
        #     alphabet[:axis],
        #     alphabet[axis:rank],
        #     alphabet[:axis],
        #     alphabet[axis:rank],
        # )

        # einsum_string_cols = '{}y{},z->{}{}zy'.format(
        #     alphabet[:axis],
        #     alphabet[axis:rank],
        #     alphabet[:axis],
        #     alphabet[axis:rank],
        # )

        # self.declare_partials(out_name, angle_name, rows=r, cols=c)
        self.declare_partials(of='*', wrt='*', method='cs')

    def compute(self, inputs, outputs):
        axis = self.options['axis']

        vec_name = self.options['vec_name']
        angle_name = self.options['angle_name']
        out_name = self.options['out_name']

        outputs[out_name] = np.concatenate([
            self.quat_real * np.cos(inputs[angle_name] / 2.),
            np.sin(inputs[angle_name] / 2.) * inputs[vec_name]
        ],
                                           axis=axis)

    # def compute_partials(self, inputs, partials):
    #     vec_name = self.options['vec_name']
    #     angle_name = self.options['angle_name']
    #     out_name = self.options['out_name']
    #     axis = self.options['axis']

    #     print(inputs[angle_name])

    #     print(self.quat_real * -0.5 * np.sin(inputs[angle_name] / 2.))

    #     comp = np.concatenate([
    #         self.quat_real * -0.5 * np.sin(inputs[angle_name] / 2.),
    #         0.5 * np.cos(inputs[angle_name] / 2.) * inputs[vec_name]
    #     ],
    #                           axis=axis)

    #     print(comp)
    #     print(comp.shape)
    #     # print(comp.reshape(24, 4, 1))

    #     partials[out_name,
    #              angle_name] = np.concatenate([
    #                  self.quat_real * -0.5 * np.sin(inputs[angle_name] / 2.),
    #                  0.5 * np.cos(inputs[angle_name] / 2.) * inputs[vec_name]
    #              ],
    #                                           axis=axis).flatten(order='F')


if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    vecshape = (3, 5)
    angle_shape = (1, 5)
    axis = 0

    prob = Problem()

    vecval = np.prod(vecshape)
    angleval = np.prod(angle_shape)
    comp = IndepVarComp()
    comp.add_output('vec', val=np.arange(vecval).reshape(vecshape))
    comp.add_output('angle', val=np.arange(angleval).reshape(angle_shape))

    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = CreateQuatComp(
        vecshape=vecshape,
        angleshape=angle_shape,
        vecval=np.arange(vecval).reshape(vecshape),
        angleval=np.arange(angleval).reshape(angle_shape),
        axis=axis,
        vec_name='vec',
        angle_name='angle',
        out_name='quat',
    )
    prob.model.add_subsystem('quat_comp', comp, promotes=['*'])

    prob.setup()
    prob.run_model()

    prob.check_partials(compact_print=True)

    prob.model.list_inputs(print_arrays=True)
    prob.model.list_outputs(print_arrays=True)