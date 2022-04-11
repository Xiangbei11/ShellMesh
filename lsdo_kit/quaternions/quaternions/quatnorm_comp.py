import numpy as np

from openmdao.api import ExplicitComponent

from omtools.utils.get_array_indices import get_array_indices

alphabet = 'abcdefghijklmnopqrstuvw'


class QuatNormComp(ExplicitComponent):
    def initialize(self):
        self.options.declare('shape', types=tuple)
        self.options.declare('val', types=np.ndarray)
        self.options.declare('axis', types=int)
        self.options.declare('in_name',
                             types=str,
                             desc="This the name of the input quaternion")
        self.options.declare('out_name',
                             types=str,
                             desc="The name of the 3x3 rotation matrix")

    def setup(self):
        shape = self.options['shape']
        val = self.options['val']
        axis = self.options['axis']
        in_name = self.options['in_name']
        out_name = self.options['out_name']

        self.add_input(in_name, shape=shape, val=val)
        self.add_output(out_name, shape=shape)

        indices = get_array_indices(*shape)

        # print(indices)

        self.shape_without_axis = shape[:axis] + shape[axis + 1:]

        ones = np.ones(4, int)

        rank = len(self.shape_without_axis)

        einsum_string_rows = '{}y{},z->{}{}yz'.format(
            alphabet[:axis],
            alphabet[axis:rank],
            alphabet[:axis],
            alphabet[axis:rank],
        )

        # print('1st:', alphabet[:axis])
        # print('2nd:', alphabet[axis:rank])

        einsum_string_cols = '{}y{},z->{}{}zy'.format(
            alphabet[:axis],
            alphabet[axis:rank],
            alphabet[:axis],
            alphabet[axis:rank],
        )

        rows = np.einsum(
            einsum_string_rows,
            indices,
            ones,
        ).flatten()

        # print(einsum_string_rows)
        # print(indices.shape)
        # print(indices)
        # print(ones.shape)
        # print(ones)

        print('rows:', rows.shape)

        cols = np.einsum(
            einsum_string_cols,
            indices,
            ones,
        ).flatten()

        # print(einsum_string_cols)
        # print(indices.shape)
        # print(indices)
        # print(ones.shape)
        # print(ones)

        print('cols:', cols.shape)

        # self.declare_partials(out_name, in_name, rows=rows, cols=cols)
        self.declare_partials(of='*', wrt='*', method='cs')

    def compute(self, inputs, outputs):
        shape = self.options['shape']
        axis = self.options['axis']
        in_name = self.options['in_name']
        out_name = self.options['out_name']

        # This computes the norm, then expands the dimension of the norm to fit the dimension of
        # the original input

        outputs[out_name] = inputs[in_name] / np.expand_dims(
            np.linalg.norm(inputs[in_name], axis=axis), axis=axis)

    # def compute_partials(self, inputs, partials):
    #     in_name = self.options['in_name']
    #     out_name = self.options['out_name']
    #     shape = self.options['shape']
    #     axis = self.options['axis']

    # partials[out_name, in_name] =


if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    shape = (2, 3, 4, 7, 5)

    prob = Problem()

    val = np.prod(shape)
    comp = IndepVarComp()
    comp.add_output('quat', val=np.random.rand(val).reshape(shape))
    prob.model.add_subsystem('ivc', comp, promotes=['*'])

    comp = QuatNormComp(
        shape=shape,
        val=np.arange(val).reshape(shape),
        axis=2,
        in_name='quat',
        out_name='quatnorm',
    )
    prob.model.add_subsystem('quatnorm_comp', comp, promotes=['*'])

    prob.setup()
    prob.run_model()
    prob.check_partials(compact_print=True)

    # prob.model.list_inputs(print_arrays=True)
    # prob.model.list_outputs(print_arrays=True)